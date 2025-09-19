// wasm-app.js - Pure C++ WebAssembly Implementation
// Simplified version with only C++ WASM support (no fallbacks)

let mcmcModule; // C++ MCMC WebAssembly module
let mcmcSamplers = []; // Array of C++ samplers for 4 chains
let mcmcAdaptiveData = []; // Array of adaptive data for 4 chains
let simulationData = null; // Store SIR simulation results
let observedData = null; // Store synthetic observed data for MCMC fitting
let customData = null; // Store user-uploaded custom data
let isUsingCustomData = false; // Track whether we're using custom data
let mcmcAlgorithm = 'adaptive'; // Track selected MCMC algorithm type
let mcmcRunning = false;
let mcmcData = {
    chains: [
        { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 },
        { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 },
        { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 },
        { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 }
    ],
    currentStep: 0,
    totalAcceptedSteps: 0,
    bestParams: null,
    bestLogPost: -Infinity
};

// SIR Model parameters
const DEFAULT_PARAMS = {
    beta: 0.7,
    gamma: 0.2,
    sigma: 0.05,
    initial_infected: 1.0,
    time_period: 200
};

// MCMC parameters
const MCMC_PARAMS = {
    adaptation_start: 50,
    adaptation_interval: 10,
    target_acceptance: 0.234,
    max_steps: 2000,
    burnin: 500
};

// Initialize the application
async function initApp() {
    try {
        // Check if Plotly is available
        if (typeof Plotly === 'undefined') {
            throw new Error('Plotly.js failed to load');
        }
        
        updateStatus('🚀 Loading C++ WebAssembly MCMC module...');
        
        // Initialize C++ MCMC module (required)
        await initializeMCMCModule();
        
        // Initialize the SIR model plot
        await setupInitialSIRPlot();
        
        // Setup event listeners
        setupEventListeners();
        
        // Hide loading screen and show main content
        document.getElementById('loading').style.display = 'none';
        document.getElementById('main-content').style.display = 'block';
        
        // Initialize MCMC button state
        updateMCMCButtonState();
        
        // Initialize intervention button state
        updateInterventionButtonState();
        
        updateStatus('✅ C++ WebAssembly MCMC ready!');
        
    } catch (error) {
        console.error('❌ Initialization failed:', error);
        document.getElementById('error-message').style.display = 'block';
        document.getElementById('error-text').textContent = `Initialization Error: ${error.message}`;
    }
}

// Initialize C++ MCMC WebAssembly module
async function initializeMCMCModule() {
    try {
        if (typeof MCMCModule === 'undefined') {
            throw new Error('C++ MCMC WebAssembly module not available. Make sure mcmc_module.js is loaded.');
        }
        
        console.log('[MCMC] Initializing C++ MCMC module...');
        mcmcModule = await MCMCModule();
        console.log('[MCMC] Module created:', mcmcModule);
        
        // Test if basic functions are available
        const testFunctions = ['create_sampler', 'create_adaptive_data', 'run_sir_simulation'];
        for (const funcName of testFunctions) {
            if (!mcmcModule.ccall || typeof mcmcModule.ccall !== 'function') {
                throw new Error('ccall function not available in WASM module');
            }
        }
        
        // Create 4 samplers for parallel chains
        mcmcSamplers = [];
        mcmcAdaptiveData = [];
        for (let i = 0; i < 4; i++) {
            const sampler = mcmcModule.ccall('create_sampler', 'number', ['number'], [12345 + i]);
            const adaptiveData = mcmcModule.ccall('create_adaptive_data', 'number', [], []);
            console.log(`[MCMC] Chain ${i}: sampler=${sampler}, adaptiveData=${adaptiveData}`);
            mcmcSamplers.push(sampler);
            mcmcAdaptiveData.push(adaptiveData);
        }
        
        // Test basic SIR simulation to verify WASM is working
        console.log('[MCMC] Testing basic SIR simulation...');
        await testSIRSimulation();
        
        console.log('[MCMC] ✅ C++ MCMC module initialized with 4 adaptive chains');
        return mcmcModule;
        
    } catch (error) {
        console.error('[MCMC] Failed to initialize C++ MCMC module:', error);
        throw error;
    }
}

// Test function to verify WASM SIR simulation works
async function testSIRSimulation() {
    try {
        const testData = await runSIRSimulation(0.7, 0.2, 1.0, 50);
        console.log('[MCMC] Test SIR simulation successful:', {
            timePoints: testData.time.length,
            maxInfected: Math.max(...testData.I).toFixed(2),
            finalSusceptible: testData.S[testData.S.length - 1].toFixed(2)
        });
    } catch (error) {
        console.error('[MCMC] Test SIR simulation failed:', error);
        throw new Error('Basic SIR simulation test failed');
    }
}

// Setup initial SIR plot
async function setupInitialSIRPlot() {
    try {
        // Run initial SIR simulation with default parameters
        const data = await runSIRSimulation(
            DEFAULT_PARAMS.beta, 
            DEFAULT_PARAMS.gamma,
            DEFAULT_PARAMS.initial_infected,
            DEFAULT_PARAMS.time_period,
            1000 // Default population size
        );
        
        simulationData = data;
        plotSIRResults(data);
        
        // Generate synthetic observed data with noise
        observedData = generateObservedData(data);
        
        console.log('[SIR] Initial plot setup complete');
        
    } catch (error) {
        console.error('[SIR] Setup error:', error);
        throw error;
    }
}

// Run SIR epidemic simulation using C++ WASM
async function runSIRSimulation(beta, gamma, initial_infected, time_period, population = 100000) {
    if (!mcmcModule || !mcmcSamplers[0]) {
        throw new Error('C++ MCMC module not initialized');
    }
    
    try {
        console.log('[SIR] Running simulation with C++ WebAssembly');
        
        // Allocate memory for output arrays
        const maxLength = time_period + 1;
        const timesPtr = mcmcModule._malloc(maxLength * 8);
        const sPtr = mcmcModule._malloc(maxLength * 8);  
        const iPtr = mcmcModule._malloc(maxLength * 8);
        const rPtr = mcmcModule._malloc(maxLength * 8);
        const lengthPtr = mcmcModule._malloc(4);
        
        try {
            // Call C++ SIR simulation
            const success = mcmcModule.ccall('run_sir_simulation', 'number', [
                'number', 'number', 'number', 'number', 'number', 'number',
                'number', 'number', 'number', 'number', 'number'
            ], [
                mcmcSamplers[0], beta, gamma, initial_infected, time_period, population,
                timesPtr, sPtr, iPtr, rPtr, lengthPtr
            ]);
            
            if (!success) {
                throw new Error('C++ SIR simulation failed');
            }
            
            // Read results
            const length = mcmcModule.getValue(lengthPtr, 'i32');
    const times = [];
    const S = [];
    const I = [];
    const R = [];
    
            for (let i = 0; i < length; i++) {
                times.push(mcmcModule.getValue(timesPtr + i * 8, 'double'));
                S.push(mcmcModule.getValue(sPtr + i * 8, 'double'));
                I.push(mcmcModule.getValue(iPtr + i * 8, 'double'));
                R.push(mcmcModule.getValue(rPtr + i * 8, 'double'));
            }
            
            console.log('[SIR] C++ WebAssembly simulation completed');
            return { time: times, S, I, R };
            
        } finally {
            // Free allocated memory
            mcmcModule._free(timesPtr);
            mcmcModule._free(sPtr);
            mcmcModule._free(iPtr);
            mcmcModule._free(rPtr);
            mcmcModule._free(lengthPtr);
        }
        
    } catch (error) {
        console.error('[SIR] C++ simulation error:', error);
        throw error;
    }
}

// Generate synthetic observed data with noise for MCMC
function generateObservedData(trueData) {
    const noise_sd = 0.05;
    const observed = {
        time: [],
        I_observed: []
    };
    
    // Sample every 10th point to simulate sparse observations
    for (let i = 0; i < trueData.time.length; i += 10) {
        if (trueData.I[i] > 0.1) { // Only observe when infection is significant
            const noise = (Math.random() - 0.5) * 2 * noise_sd * trueData.I[i];
            observed.time.push(trueData.time[i]);
            observed.I_observed.push(Math.max(0.01, trueData.I[i] + noise));
        }
    }
    
    console.log(`[MCMC] Generated ${observed.time.length} synthetic observations`);
    return observed;
}

// Run simulation button handler
async function runSimulation() {
    try {
        updateStatus('Running epidemic simulation...');
        
        const beta = parseFloat(document.getElementById('beta').value);
        const gamma = parseFloat(document.getElementById('gamma').value);
        const initial_infected = parseFloat(document.getElementById('initial_infected').value);
        const time_period = parseInt(document.getElementById('time_period').value);
        const population_size = parseInt(document.getElementById('population_size').value);
        
        // Run the SIR simulation using C++ WASM
        simulationData = await runSIRSimulation(beta, gamma, initial_infected, time_period, population_size);
        
        // Plot the results
        plotSIRResults(simulationData);
        
        // Generate synthetic observed data with noise
        observedData = generateObservedData(simulationData);
        
        // Enable MCMC button
        document.getElementById('mcmc-btn').disabled = false;
        updateMCMCButtonState();
        
        updateStatus('Data generated - ready for MCMC fitting');
        
    } catch (error) {
        console.error('[SIR] Simulation error:', error);
        updateStatus('Simulation failed');
    }
}

// Plot SIR simulation results
function plotSIRResults(data, animateToFrame = null) {
    const endFrame = animateToFrame || data.time.length - 1;
    
    const susceptibleTrace = {
        x: data.time.slice(0, endFrame + 1),
        y: data.S.slice(0, endFrame + 1),
        type: 'scatter',
        mode: 'lines',
        line: { color: '#3498db', width: 3 },
        name: 'Susceptible'
    };
    
    const infectedTrace = {
        x: data.time.slice(0, endFrame + 1),
        y: data.I.slice(0, endFrame + 1),
        type: 'scatter',
        mode: 'lines',
        line: { color: '#e74c3c', width: 3 },
        name: 'Infected'
    };
    
    const recoveredTrace = {
        x: data.time.slice(0, endFrame + 1),
        y: data.R.slice(0, endFrame + 1),
        type: 'scatter',
        mode: 'lines',
        line: { color: '#27ae60', width: 3 },
        name: 'Recovered'
    };
    
    // Only show the basic SIR simulation data
    const traces = [susceptibleTrace, infectedTrace, recoveredTrace];
    
    const layout = {
        title: 'SIR Epidemic Model',
        xaxis: { title: 'Time (days)' },
        yaxis: { title: 'Percentage of Population' },
        showlegend: true,
        height: 400
    };
    
    Plotly.newPlot('plotContainer', traces, layout, { responsive: true, displayModeBar: false });
}

// Start MCMC sampling with C++ WebAssembly
async function startMCMC() {
    if (mcmcRunning) return;
    
    // Check if we have data to run MCMC on
    if (!observedData || !observedData.time || observedData.time.length === 0) {
        showDataWarning();
        return;
    }
    
    mcmcRunning = true;
    document.getElementById('mcmc-btn').disabled = true;
    document.getElementById('stop-mcmc-btn').disabled = false;
    
    // Reset MCMC data for 4 chains
    mcmcData = {
        chains: [
            { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 },
            { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 },
            { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 },
            { traces: { beta: [], gamma: [], sigma: [] }, logLikelihood: [], currentStep: 0, acceptedSteps: 0 }
        ],
        currentStep: 0,
        totalAcceptedSteps: 0,
        bestParams: null,
        bestLogPost: -Infinity,
        chainParams: [
            { beta: 0.5, gamma: 0.15, sigma: 0.03 },
            { beta: 0.6, gamma: 0.18, sigma: 0.04 },
            { beta: 0.8, gamma: 0.22, sigma: 0.06 },
            { beta: 0.7, gamma: 0.20, sigma: 0.05 }
        ],
        chainLogPost: [-Infinity, -Infinity, -Infinity, -Infinity]
    };
    
    // Reset adaptive data
    for (let i = 0; i < 4; i++) {
        mcmcModule.ccall('reset_adaptive_data', 'void', ['number'], [mcmcAdaptiveData[i]]);
    }
    
    document.getElementById('mcmc-progress').textContent = 'Step 0 / 0 (0%)';
    document.getElementById('progress-bar').style.width = '0%';
    document.getElementById('current-loglik').textContent = '-';
    document.getElementById('accept-rate').textContent = '-';
    document.getElementById('current-params').textContent = 'β:- γ:- σ:-';
    
    // Initialize convergence status
    const indicatorElement = document.getElementById('convergence-indicator');
    if (indicatorElement) {
        indicatorElement.innerHTML = '⏳ <span style="color: #6c757d;">Running...</span>';
    }
    
    updateStatus('Starting MCMC sampling with C++ WebAssembly...');
    
    try {
        // Read values from sliders
        const maxSteps = parseInt(document.getElementById('mcmc_steps').value);
        const burnin = parseInt(document.getElementById('burnin').value);
        
        console.log(`[MCMC] Starting ${maxSteps} steps with ${burnin} burnin (4 chains) using ${mcmcAlgorithm} algorithm`);
        await runMCMCLoop(mcmcData.chainParams, maxSteps, burnin);
        
    } catch (error) {
        console.error('[MCMC] Error:', error);
        updateStatus('MCMC sampling failed');
        stopMCMC();
    }
}

// MCMC sampling loop for 4 adaptive chains using C++ WASM
async function runMCMCLoop(initialParams, maxSteps, burnin) {
    for (let step = 0; step < maxSteps && mcmcRunning; step++) {
        mcmcData.currentStep = step;
        
        const chainPromises = [];
        
        // Run all 4 chains in parallel using C++ WASM
        for (let chainId = 0; chainId < 4; chainId++) {
            chainPromises.push(
                mcmcStepCPPAdaptive(chainId, mcmcData.chainParams[chainId], mcmcData.chainLogPost[chainId])
            );
        }
        
        try {
            const results = await Promise.all(chainPromises);
        let stepBestParams = null;
        
            // Process results from all chains
        for (let chainId = 0; chainId < 4; chainId++) {
                const result = results[chainId];
                
                // Update chain data
            mcmcData.chainParams[chainId] = result.params;
                mcmcData.chainLogPost[chainId] = result.log_posterior;
                
                // Store traces
                mcmcData.chains[chainId].traces.beta.push(result.params.beta);
                mcmcData.chains[chainId].traces.gamma.push(result.params.gamma);
                mcmcData.chains[chainId].traces.sigma.push(result.params.sigma);
                mcmcData.chains[chainId].logLikelihood.push(result.log_posterior);
                mcmcData.chains[chainId].currentStep = step;
                
                if (result.accepted) {
                mcmcData.chains[chainId].acceptedSteps++;
                    mcmcData.totalAcceptedSteps++;
                }
                
                // Track best parameters
                if (result.log_posterior > mcmcData.bestLogPost) {
                    mcmcData.bestLogPost = result.log_posterior;
                    mcmcData.bestParams = { ...result.params };
                    stepBestParams = { ...result.params };
                }
                
                // Limit trace length for memory
                const MAX_TRACE_LENGTH = 20000;
            if (mcmcData.chains[chainId].traces.beta.length > MAX_TRACE_LENGTH) {
                const excess = mcmcData.chains[chainId].traces.beta.length - MAX_TRACE_LENGTH;
                mcmcData.chains[chainId].traces.beta = mcmcData.chains[chainId].traces.beta.slice(excess);
                mcmcData.chains[chainId].traces.gamma = mcmcData.chains[chainId].traces.gamma.slice(excess);
                mcmcData.chains[chainId].traces.sigma = mcmcData.chains[chainId].traces.sigma.slice(excess);
                mcmcData.chains[chainId].logLikelihood = mcmcData.chains[chainId].logLikelihood.slice(excess);
                    console.log(`[MCMC] Trimmed traces for chain ${chainId} (keeping last ${MAX_TRACE_LENGTH})`);
                }
                
                // Enable adaptive features only if using adaptive algorithm
                if (mcmcAlgorithm === 'adaptive') {
                    if (step === MCMC_PARAMS.adaptation_start) {
                        const success = mcmcModule.ccall('enable_adaptive_sampling', 'number', ['number'], [mcmcAdaptiveData[chainId]]);
                        if (success) {
                            console.log(`[MCMC] Chain ${chainId} enabled adaptive sampling at step ${step}`);
                        }
                    } else if (step > MCMC_PARAMS.adaptation_start && step % MCMC_PARAMS.adaptation_interval === 0) {
                        // Update empirical statistics periodically
                        mcmcModule.ccall('update_empirical_stats', 'number', ['number', 'number'], [mcmcSamplers[chainId], mcmcAdaptiveData[chainId]]);
                    }
                }
            }
            
            // Update UI every 10 steps
            if (step % 10 === 0 || step < 100) {
                updateMCMCDisplays(step, maxSteps);
                
                // Update MCMC plots every 50 steps (less frequent to avoid performance issues)
                if (step % 50 === 0 && step > 0) {
                    await updateMCMCPlots();
                }
                
                if (stepBestParams && simulationData) {
                    try {
                        // Show prediction with best parameters
                        console.log(`[MCMC] Plotting step ${step} with params: β=${stepBestParams.beta.toFixed(3)}, γ=${stepBestParams.gamma.toFixed(3)}, σ=${stepBestParams.sigma.toFixed(3)}`);
                        
                        const predictedData = await runSIRSimulation(
                            stepBestParams.beta,
                            stepBestParams.gamma,
                            parseFloat(document.getElementById('initial_infected').value),
                            parseInt(document.getElementById('time_period').value),
                            parseInt(document.getElementById('population_size').value)
                        );
                        
                        if (predictedData && predictedData.time && predictedData.time.length > 0) {
                            plotMCMCResults(simulationData, predictedData);
                            console.log(`[MCMC] Successfully plotted step ${step}`);
                } else {
                            console.warn(`[MCMC] Invalid predicted data at step ${step}`);
                        }
                    } catch (plotError) {
                        console.error(`[MCMC] Plotting error at step ${step}:`, plotError);
                    }
                }
            }
            
        } catch (error) {
            console.error(`[MCMC] Error in step ${step}:`, error);
        }
        
        // Small delay for UI responsiveness
        if (step % 50 === 0) {
            await new Promise(resolve => setTimeout(resolve, 1));
        }
    }
    
    if (mcmcRunning) {
        console.log('[MCMC] ✅ MCMC sampling completed');
        updateStatus('MCMC sampling completed');
    mcmcRunning = false;
    document.getElementById('mcmc-btn').disabled = false;
    document.getElementById('stop-mcmc-btn').disabled = true;
    
    // Final plot update
    await updateMCMCPlots();
    
    // Final convergence check
    updateConvergenceDiagnostics();
    
    // Update download button state
    updateDownloadButtonState();
    
    // Update intervention button state
    updateInterventionButtonState();
    }
}

// C++ MCMC step function (handles both standard and adaptive)
async function mcmcStepCPPAdaptive(chainId, currentParams, currentLogPost) {
    if (!mcmcModule || !mcmcSamplers[chainId] || !mcmcAdaptiveData[chainId]) {
        throw new Error('C++ MCMC module not initialized');
    }
    
    try {
        // Allocate memory for observed data
        const nObs = observedData.time.length;
        const obsTimesPtr = mcmcModule._malloc(nObs * 8);
        const obsInfectedPtr = mcmcModule._malloc(nObs * 8);
        
        // Allocate memory for output parameters
        const newBetaPtr = mcmcModule._malloc(8);
        const newGammaPtr = mcmcModule._malloc(8);
        const newSigmaPtr = mcmcModule._malloc(8);
        const newLogPostPtr = mcmcModule._malloc(8);
        const acceptedPtr = mcmcModule._malloc(4);
        const usedAdaptivePtr = mcmcModule._malloc(4);
        const acceptanceRatePtr = mcmcModule._malloc(8);
        const scalingFactorPtr = mcmcModule._malloc(8);
        
        // Allocate memory for predicted data
        const maxPredLength = parseInt(document.getElementById('time_period').value) + 1;
        const predTimesPtr = mcmcModule._malloc(maxPredLength * 8);
        const predInfectedPtr = mcmcModule._malloc(maxPredLength * 8);
        const predLengthPtr = mcmcModule._malloc(4);
        
        try {
            // Copy observed data to WASM memory
            for (let i = 0; i < nObs; i++) {
            mcmcModule.setValue(obsTimesPtr + i * 8, observedData.time[i], 'double');
            mcmcModule.setValue(obsInfectedPtr + i * 8, observedData.I_observed[i], 'double');
        }
        
            // Get form values
            const initial_infected_pct = parseFloat(document.getElementById('initial_infected').value);
        const time_period = parseInt(document.getElementById('time_period').value);
        
            console.log(`[MCMC] Chain ${chainId} using ${mcmcAlgorithm} algorithm`);
            let success;
            
            if (mcmcAlgorithm === 'adaptive') {
                // Call adaptive MCMC step function
                success = mcmcModule.ccall('mcmc_step_adaptive', 'number', [
                    'number', 'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number'
                ], [
                    mcmcSamplers[chainId], mcmcAdaptiveData[chainId],
                    currentParams.beta, currentParams.gamma, currentParams.sigma,
                    currentLogPost,
                    obsTimesPtr, obsInfectedPtr, nObs,
                    0.7, 0.3,  // beta prior mean, std
                    0.2, 0.1,  // gamma prior mean, std  
                    2.0, 1.0,  // sigma prior a, b
                    initial_infected_pct, time_period,
                    newBetaPtr, newGammaPtr, newSigmaPtr, newLogPostPtr, acceptedPtr,
                    usedAdaptivePtr, acceptanceRatePtr, scalingFactorPtr,
                    predTimesPtr, predInfectedPtr, predLengthPtr
                ]);
            } else {
                // Call standard MCMC step function with fixed proposal standard deviations
                success = mcmcModule.ccall('mcmc_step', 'number', [
                    'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number', 
                    'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number', 'number', 'number',
                    'number', 'number', 'number'
                ], [
                    mcmcSamplers[chainId],
                    currentParams.beta, currentParams.gamma, currentParams.sigma,
                    currentLogPost,
                    obsTimesPtr, obsInfectedPtr, nObs,
                    0.7, 0.3,  // beta prior mean, std
                    0.2, 0.1,  // gamma prior mean, std  
                    2.0, 1.0,  // sigma prior a, b
                    0.02, 0.01, 0.005,  // beta, gamma, sigma proposal standard deviations
                    initial_infected_pct, time_period,
                    newBetaPtr, newGammaPtr, newSigmaPtr, newLogPostPtr, acceptedPtr,
                    predTimesPtr, predInfectedPtr, predLengthPtr
                ]);
            }
        
        if (!success) {
                throw new Error(`C++ ${mcmcAlgorithm} MCMC step failed`);
        }
        
            // Read results based on algorithm
            const result = {
                params: {
                    beta: mcmcModule.getValue(newBetaPtr, 'double'),
                    gamma: mcmcModule.getValue(newGammaPtr, 'double'),
                    sigma: mcmcModule.getValue(newSigmaPtr, 'double')
                },
                log_posterior: mcmcModule.getValue(newLogPostPtr, 'double'),
                accepted: mcmcModule.getValue(acceptedPtr, 'i32') === 1,
                used_adaptive: mcmcAlgorithm === 'adaptive' ? mcmcModule.getValue(usedAdaptivePtr, 'i32') === 1 : false,
                acceptance_rate: mcmcAlgorithm === 'adaptive' ? mcmcModule.getValue(acceptanceRatePtr, 'double') : 0.0,
                scaling_factor: mcmcAlgorithm === 'adaptive' ? mcmcModule.getValue(scalingFactorPtr, 'double') : 1.0
            };
            
            return result;
            
        } finally {
            // Free all allocated memory
            mcmcModule._free(obsTimesPtr);
            mcmcModule._free(obsInfectedPtr);
            mcmcModule._free(newBetaPtr);
            mcmcModule._free(newGammaPtr);
            mcmcModule._free(newSigmaPtr);
            mcmcModule._free(newLogPostPtr);
            mcmcModule._free(acceptedPtr);
            mcmcModule._free(usedAdaptivePtr);
            mcmcModule._free(acceptanceRatePtr);
            mcmcModule._free(scalingFactorPtr);
            mcmcModule._free(predTimesPtr);
            mcmcModule._free(predInfectedPtr);
            mcmcModule._free(predLengthPtr);
        }
        
    } catch (error) {
        console.error(`[MCMC] C++ adaptive step failed for chain ${chainId}:`, error);
        throw error;
    }
}

// Update MCMC displays
function updateMCMCDisplays(step, maxSteps) {
    const progress = ((step + 1) / maxSteps * 100).toFixed(1);
    document.getElementById('mcmc-progress').textContent = `Step ${step + 1} / ${maxSteps} (${progress}%)`;
    document.getElementById('progress-bar').style.width = progress + '%';
    
    // Calculate overall acceptance rate
    let totalSteps = 0;
        let totalAccepted = 0;
    for (let i = 0; i < 4; i++) {
        totalSteps += mcmcData.chains[i].currentStep + 1;
        totalAccepted += mcmcData.chains[i].acceptedSteps;
    }
    const acceptRate = totalSteps > 0 ? (totalAccepted / totalSteps * 100).toFixed(1) : '0.0';
    document.getElementById('accept-rate').textContent = acceptRate + '%';
    
            if (mcmcData.bestParams) {
        const paramsText = `β:${mcmcData.bestParams.beta.toFixed(2)} γ:${mcmcData.bestParams.gamma.toFixed(2)} σ:${mcmcData.bestParams.sigma.toFixed(2)}`;
        document.getElementById('current-params').textContent = paramsText;
        document.getElementById('current-loglik').textContent = mcmcData.bestLogPost.toFixed(2);
    }
    
    // Show algorithm status
    const adaptiveElement = document.getElementById('cov-status');
    if (adaptiveElement) {
        if (mcmcAlgorithm === 'adaptive') {
            if (step > MCMC_PARAMS.adaptation_start) {
                adaptiveElement.innerHTML = `Active: 4/4 chains<br>Adaptive covariance + Robbins-Monro scaling`;
            } else {
                const stepsUntil = MCMC_PARAMS.adaptation_start - step;
                adaptiveElement.innerHTML = `Inactive<br>${stepsUntil} steps until adaptive`;
            }
        } else {
            adaptiveElement.innerHTML = `Standard M-H<br>Fixed proposal covariance`;
        }
    }
}

// Plot MCMC results showing original vs fitted
function plotMCMCResults(originalData, fittedData) {
    console.log('[MCMC] plotMCMCResults called with:', {
        isCustomData: isUsingCustomData,
        originalDataLength: originalData?.time?.length || 0,
        fittedDataLength: fittedData?.time?.length || 0,
        originalSample: originalData?.I?.slice(0, 3),
        fittedSample: fittedData?.I?.slice(0, 3)
    });
    
    let traces = [];
    
    if (isUsingCustomData && customData) {
        // For custom data, show the custom data and model fit (if available)
        traces.push({
            x: customData.time,
            y: customData.incidence,
            type: 'scatter',
            mode: 'lines+markers',
            line: { color: '#e74c3c', width: 3 },
            marker: { color: '#e74c3c', size: 6 },
            name: 'Custom Data'
        });
        
        // If we have a fitted model for the custom data time points, show it
        if (fittedData && fittedData.time && fittedData.I) {
            traces.push({
                x: fittedData.time,
                y: fittedData.I,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#9b59b6', width: 2, dash: 'dash' },
                name: 'MCMC Model Fit'
            });
        }
        
        // Add observed data points (which are the same as custom data)
        if (observedData) {
            traces.push({
                x: observedData.time,
                y: observedData.I_observed,
                type: 'scatter',
                mode: 'markers',
                marker: { color: '#f39c12', size: 8, symbol: 'circle-open' },
                name: 'Data Points'
            });
        }
        
    } else {
        // Original behavior for simulated data
        if (!originalData || !fittedData || !originalData.time || !fittedData.time) {
            console.error('[MCMC] Invalid data passed to plotMCMCResults');
            return;
        }
        
        traces = [
            {
                x: originalData.time,
                y: originalData.S,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#3498db', width: 2, dash: 'dash' },
                name: 'True Susceptible'
            },
            {
                x: originalData.time,
                y: originalData.I,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#e74c3c', width: 2, dash: 'dash' },
                name: 'True Infected'
            },
            {
                x: originalData.time,
                y: originalData.R,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#27ae60', width: 2, dash: 'dash' },
                name: 'True Recovered'
            },
            {
                x: fittedData.time,
                y: fittedData.S,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#3498db', width: 3 },
                name: 'Fitted Susceptible'
            },
            {
                x: fittedData.time,
                y: fittedData.I,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#e74c3c', width: 3 },
                name: 'Fitted Infected'
            },
            {
                x: fittedData.time,
                y: fittedData.R,
                type: 'scatter',
                mode: 'lines',
                line: { color: '#27ae60', width: 3 },
                name: 'Fitted Recovered'
            }
        ];
        
        // Add observed data points
        if (observedData) {
            traces.push({
                x: observedData.time,
                y: observedData.I_observed,
                type: 'scatter',
                mode: 'markers',
                marker: { color: '#f39c12', size: 8, symbol: 'circle-open' },
                name: 'Observed Data'
            });
        }
    }
    
    const layout = {
        title: isUsingCustomData ? 'MCMC Fitting - Custom Data' : 'MCMC Bayesian Fitting - C++ WebAssembly',
        xaxis: { title: isUsingCustomData ? 'Time' : 'Time (days)' },
        yaxis: { title: isUsingCustomData ? 'Incidence' : 'Percentage of Population' },
        showlegend: true,
        height: 400
    };
    
    try {
        console.log('[MCMC] Calling Plotly.newPlot with', traces.length, 'traces');
        Plotly.newPlot('plotContainer', traces, layout, { responsive: true, displayModeBar: false });
        console.log('[MCMC] ✅ Plot updated successfully');
    } catch (plotError) {
        console.error('[MCMC] Plotly error:', plotError);
        throw plotError;
    }
}

// Stop MCMC sampling
function stopMCMC() {
    mcmcRunning = false;
    document.getElementById('mcmc-btn').disabled = false;
    document.getElementById('stop-mcmc-btn').disabled = true;
    updateStatus('MCMC sampling stopped');
    console.log('[MCMC] Sampling stopped by user');
}

// Update R0 display (removed from UI)
function updateR0Display() {
    // R0 display removed from interface - function kept for compatibility
    return;
}

// Update MCMC button state based on data availability
function updateMCMCButtonState() {
    const mcmcBtn = document.getElementById('mcmc-btn');
    const progressElement = document.getElementById('mcmc-progress');
    
    if (!mcmcBtn) return;
    
    if (!observedData || !observedData.time || observedData.time.length === 0) {
        // No data available
        mcmcBtn.disabled = true;
        mcmcBtn.textContent = 'No Data - Simulate or Upload First';
        mcmcBtn.style.backgroundColor = '#e74c3c';
        mcmcBtn.style.color = 'white';
        mcmcBtn.style.borderColor = '#e74c3c';
        
        if (progressElement) {
            progressElement.innerHTML = '<span style="color: #e74c3c; font-weight: bold;">⚠️ No Data Available</span>';
        }
    } else {
        // Data available
        mcmcBtn.disabled = false;
        mcmcBtn.textContent = 'Start';
        mcmcBtn.style.backgroundColor = '';
        mcmcBtn.style.color = '';
        mcmcBtn.style.borderColor = '';
        
        if (progressElement) {
            progressElement.textContent = 'Step 0 / 0 (0%)';
        }
    }
}

// Show data warning when trying to start MCMC without data
function showDataWarning() {
    const warningMessage = '⚠️ No data available for MCMC fitting. Please either:\n• Click "Simulate Data" to generate synthetic data, or\n• Upload custom data using the file upload button';
    
    // Show warning in status
    updateStatus(warningMessage);
    
    // Show visual warning in the MCMC progress area
    const progressElement = document.getElementById('mcmc-progress');
    if (progressElement) {
        progressElement.innerHTML = '<span style="color: #e74c3c; font-weight: bold;">⚠️ No Data Available</span>';
    }
    
    // Temporarily disable the MCMC button and show why
    const mcmcBtn = document.getElementById('mcmc-btn');
    if (mcmcBtn) {
        mcmcBtn.disabled = true;
        mcmcBtn.textContent = 'No Data - Simulate or Upload First';
        mcmcBtn.style.backgroundColor = '#e74c3c';
        mcmcBtn.style.color = 'white';
        mcmcBtn.style.borderColor = '#e74c3c';
    }
    
    // Re-enable button after 3 seconds
    setTimeout(() => {
        if (mcmcBtn) {
            mcmcBtn.disabled = false;
            mcmcBtn.textContent = 'Start';
            mcmcBtn.style.backgroundColor = '';
            mcmcBtn.style.color = '';
            mcmcBtn.style.borderColor = '';
        }
        if (progressElement) {
            progressElement.textContent = 'Step 0 / 0 (0%)';
        }
    }, 3000);
    
    console.log('[MCMC] Warning: No data available for MCMC fitting');
}

// Update status message
function updateStatus(message) {
    const statusElement = document.getElementById('status');
    const mathBackendElement = document.getElementById('math-backend');
    
    if (statusElement) {
        statusElement.textContent = message;
    }
    
    // Update math backend display
    if (mathBackendElement) {
        mathBackendElement.textContent = 'C++ WebAssembly (High Performance)';
    }
    
    console.log(message);
}

// Reset model to initial state
function resetModel() {
    // Reset to default parameters
    document.getElementById('beta').value = DEFAULT_PARAMS.beta;
    document.getElementById('gamma').value = DEFAULT_PARAMS.gamma;
    document.getElementById('initial_infected').value = DEFAULT_PARAMS.initial_infected;
    document.getElementById('time_period').value = DEFAULT_PARAMS.time_period;
    document.getElementById('population_size').value = 1000;
    
    // Update displays
    document.getElementById('beta-value').textContent = DEFAULT_PARAMS.beta.toFixed(2);
    document.getElementById('gamma-value').textContent = DEFAULT_PARAMS.gamma.toFixed(2);
    document.getElementById('initial_infected-value').textContent = DEFAULT_PARAMS.initial_infected.toFixed(1);
    document.getElementById('time_period-value').textContent = DEFAULT_PARAMS.time_period.toString();
    document.getElementById('population_size-value').textContent = '1000';
    
    // Reset custom data state
    customData = null;
    isUsingCustomData = false;
    updateUploadStatus('');
    
    // Clear file input
    const fileInput = document.getElementById('data-file-input');
    if (fileInput) {
        fileInput.value = '';
    }
    
    // Run simulation with default parameters
    runSimulation();
    
    // Update MCMC button state after reset
    updateMCMCButtonState();
    
    console.log('[UI] Model reset to defaults (including custom data)');
}

// Setup event listeners
function setupEventListeners() {
    document.getElementById('simulate-btn').addEventListener('click', runSimulation);
    document.getElementById('mcmc-btn').addEventListener('click', startMCMC);
    document.getElementById('stop-mcmc-btn').addEventListener('click', stopMCMC);
    document.getElementById('reset-btn').addEventListener('click', resetModel);
    
    // Parameter change listeners with value display updates
    const parameterSliders = [
        { id: 'beta', decimals: 2 },
        { id: 'gamma', decimals: 2 },
        { id: 'initial_infected', decimals: 1 },
        { id: 'time_period', decimals: 0 },
        { id: 'population_size', decimals: 0 },
        { id: 'mcmc_steps', decimals: 0 },
        { id: 'burnin', decimals: 0 }
    ];
    
    parameterSliders.forEach(param => {
        const slider = document.getElementById(param.id);
        const valueDisplay = document.getElementById(param.id + '-value');
        
        if (slider && valueDisplay) {
            slider.addEventListener('input', function() {
                const value = parseFloat(this.value);
                valueDisplay.textContent = value.toFixed(param.decimals);
                
                // R0 display removed from interface
                
                console.log(`[UI] ${param.id} updated to ${value}`);
            });
            
            // Initialize value display
            const initialValue = parseFloat(slider.value);
            valueDisplay.textContent = initialValue.toFixed(param.decimals);
        } else {
            console.warn(`[UI] Missing slider or value display for ${param.id}`);
        }
    });
    
    // R0 display removed from interface
    
    // Add event listener for include burn-in checkbox
    const includeBurninCheckbox = document.getElementById('include-burnin');
    if (includeBurninCheckbox) {
        includeBurninCheckbox.addEventListener('change', async function() {
            console.log(`[UI] Include burn-in changed to ${this.checked}`);
            await updateMCMCPlots(); // Update plots when burn-in setting changes
        });
    }
    
    // Add event listener for algorithm selection
    const algorithmSelect = document.getElementById('mcmc_algorithm');
    if (algorithmSelect) {
        console.log('[UI] Setting up algorithm dropdown event listener');
        algorithmSelect.addEventListener('change', function() {
            mcmcAlgorithm = this.value;
            console.log(`[UI] MCMC algorithm changed to ${mcmcAlgorithm}`);
            updateAlgorithmInfo(mcmcAlgorithm);
            updateAlgorithmStatusDisplay(mcmcAlgorithm);
        });
        
        // Initialize algorithm info with current value
        mcmcAlgorithm = algorithmSelect.value; // Make sure we get the actual dropdown value
        updateAlgorithmInfo(mcmcAlgorithm);
        updateAlgorithmStatusDisplay(mcmcAlgorithm);
        console.log(`[UI] Algorithm dropdown initialized with: ${mcmcAlgorithm}`);
    } else {
        console.warn('[UI] Algorithm dropdown element not found!');
    }
    
    // Add event listener for thinning factor slider
    const thinningSlider = document.getElementById('thinning-factor');
    const thinningValueDisplay = document.getElementById('thinning-factor-value');
    if (thinningSlider && thinningValueDisplay) {
        thinningSlider.addEventListener('input', function() {
            const value = parseInt(this.value);
            thinningValueDisplay.textContent = value;
            console.log(`[UI] Thinning factor updated to ${value}`);
        });
        
        // Initialize value display
        const initialValue = parseInt(thinningSlider.value);
        thinningValueDisplay.textContent = initialValue;
    }
    
    // Add event listeners for intervention controls
    const interventionSliders = [
        { id: 'vaccine_coverage', suffix: '%' },
        { id: 'vaccine_start', suffix: '' },
        { id: 'vaccine_end', suffix: '' }
    ];
    
    interventionSliders.forEach(param => {
        const slider = document.getElementById(param.id);
        const valueDisplay = document.getElementById(param.id + '-value');
        
        if (slider && valueDisplay) {
            slider.addEventListener('input', function() {
                const value = parseFloat(this.value);
                valueDisplay.textContent = value + param.suffix;
                console.log(`[UI] ${param.id} updated to ${value}`);
                
                // Update uptake plot for coverage, start, and end time changes
                if (['vaccine_coverage', 'vaccine_start', 'vaccine_end'].includes(param.id)) {
                    plotVaccinationUptake();
                }
            });
            
            // Initialize value display
            const initialValue = parseFloat(slider.value);
            valueDisplay.textContent = initialValue + param.suffix;
        }
    });
    
    // Add event listeners for waning controls
    const waningSliders = [
        { id: 'waning_protection', suffix: '%', decimals: 0 },
        { id: 'waning_rate', suffix: '', decimals: 3 }
    ];
    
    waningSliders.forEach(param => {
        const slider = document.getElementById(param.id);
        const valueDisplay = document.getElementById(param.id + '-value');
        
        if (slider && valueDisplay) {
            slider.addEventListener('input', function() {
                const value = parseFloat(this.value);
                // Special handling for waning rate to show 0.000 when zero
                if (param.id === 'waning_rate' && value === 0) {
                    valueDisplay.textContent = '0.000';
                } else {
                    valueDisplay.textContent = value.toFixed(param.decimals) + param.suffix;
                }
                console.log(`[UI] ${param.id} updated to ${value}`);
                
                // Update waning plot
                plotWaningCurve();
            });
            
            // Initialize value display
            const initialValue = parseFloat(slider.value);
            if (param.id === 'waning_rate' && initialValue === 0) {
                valueDisplay.textContent = '0.000';
            } else {
                valueDisplay.textContent = initialValue.toFixed(param.decimals) + param.suffix;
            }
        }
    });
    
    // Add event listener for waning distribution selection
    const waningDistributionSelect = document.getElementById('waning_distribution');
    if (waningDistributionSelect) {
        waningDistributionSelect.addEventListener('change', function() {
            console.log(`[UI] Waning distribution changed to ${this.value}`);
            updateWaningDistributionInfo();
            plotWaningCurve();
        });
        
        // Initialize distribution info
        updateWaningDistributionInfo();
    }
    
    // Initialize plots
    plotVaccinationUptake();
    plotWaningCurve();
    
    console.log('[UI] ✅ Event listeners setup complete');
}

// Plot MCMC trace plots
function plotTraces() {
    if (!mcmcData || mcmcData.chains.length === 0) return;
    
    const traces = [];
    const includeBurnin = document.getElementById('include-burnin')?.checked !== false;
    const burnin = includeBurnin ? 0 : parseInt(document.getElementById('burnin').value || 500);
    
    // Create traces for each parameter across all chains
    ['beta', 'gamma', 'sigma'].forEach((param, paramIndex) => {
        mcmcData.chains.forEach((chain, chainIndex) => {
            if (chain.traces[param] && chain.traces[param].length > burnin) {
                const data = chain.traces[param].slice(burnin);
                traces.push({
                    x: Array.from({length: data.length}, (_, i) => i + burnin),
                    y: data,
                    type: 'scatter',
                    mode: 'lines',
                    name: `${param} Chain ${chainIndex + 1}`,
                    line: { width: 1 },
                    yaxis: paramIndex === 0 ? 'y' : paramIndex === 1 ? 'y2' : 'y3'
                });
            }
        });
    });
    
    const layout = {
        title: 'Parameter Traces',
        xaxis: { title: 'MCMC Step' },
        yaxis: { title: 'β (Transmission Rate)', domain: [0.7, 1] },
        yaxis2: { title: 'γ (Recovery Rate)', domain: [0.35, 0.65] },
        yaxis3: { title: 'σ (Noise)', domain: [0, 0.3] },
        height: 250,
        margin: { l: 50, r: 30, t: 30, b: 30 },
        showlegend: false
    };
    
    Plotly.newPlot('tracePlot', traces, layout, { responsive: true, displayModeBar: false });
}

// Plot posterior density plots (faceted)
function plotPosteriors() {
    if (!mcmcData || mcmcData.chains.length === 0) return;
    
    const traces = [];
    const includeBurnin = document.getElementById('include-burnin')?.checked !== false;
    const burnin = includeBurnin ? 0 : parseInt(document.getElementById('burnin').value || 500);
    
    // Create separate faceted plots for each parameter
    ['beta', 'gamma', 'sigma'].forEach((param, paramIndex) => {
        const allData = [];
        mcmcData.chains.forEach(chain => {
            if (chain.traces[param] && chain.traces[param].length > burnin) {
                allData.push(...chain.traces[param].slice(burnin));
            }
        });
        
        if (allData.length > 0) {
            traces.push({
                x: allData,
                type: 'histogram',
                name: param,
                opacity: 0.8,
                histnorm: 'probability density',
                nbinsx: 25,
                yaxis: paramIndex === 0 ? 'y' : paramIndex === 1 ? 'y2' : 'y3',
                xaxis: paramIndex === 0 ? 'x' : paramIndex === 1 ? 'x2' : 'x3',
                marker: {
                    color: paramIndex === 0 ? '#3498db' : paramIndex === 1 ? '#e74c3c' : '#27ae60'
                }
            });
        }
    });
    
    const layout = {
        title: 'Posterior Densities',
        xaxis: { title: 'β (Transmission Rate)', domain: [0, 0.32] },
        xaxis2: { title: 'γ (Recovery Rate)', domain: [0.34, 0.66] },
        xaxis3: { title: 'σ (Noise)', domain: [0.68, 1] },
        yaxis: { title: 'Density', domain: [0, 1] },
        yaxis2: { title: '', domain: [0, 1], anchor: 'x2' },
        yaxis3: { title: '', domain: [0, 1], anchor: 'x3' },
        height: 250,
        margin: { l: 50, r: 30, t: 30, b: 30 },
        showlegend: false
    };
    
    Plotly.newPlot('posteriorPlot', traces, layout, { responsive: true, displayModeBar: false });
}

// Plot predictive fit
async function plotPredictive() {
    if (!mcmcData.bestParams || !observedData) return;
    
    const traces = [];
    
    try {
        if (isUsingCustomData && customData) {
            // For custom data, we don't have full SIR simulation - just show fitted line if possible
            // The MCMC will have fitted parameters to match the incidence data directly
            traces.push({
                x: customData.time,
                y: customData.incidence,
                type: 'scatter',
                mode: 'lines',
                name: 'Custom Data',
                line: { color: '#e74c3c', width: 2 }
            });
        } else if (simulationData) {
            // Generate prediction with best parameters for simulated data
            const predictedData = await runSIRSimulation(
                mcmcData.bestParams.beta,
                mcmcData.bestParams.gamma,
                parseFloat(document.getElementById('initial_infected').value),
                parseInt(document.getElementById('time_period').value),
                parseInt(document.getElementById('population_size').value)
            );
            
            // Add fitted infected curve
            if (predictedData && predictedData.time) {
                traces.push({
                    x: predictedData.time,
                    y: predictedData.I,
                    type: 'scatter',
                    mode: 'lines',
                    name: 'Fitted',
                    line: { color: '#e74c3c', width: 2 }
                });
            }
        }
    } catch (error) {
        console.warn('[MCMC] Could not generate prediction for predictive plot:', error);
    }
    
    // Add observed data points (either custom or synthetic)
    const dataName = isUsingCustomData ? 'Custom Points' : 'Observed';
    traces.push({
        x: observedData.time,
        y: observedData.I_observed,
        type: 'scatter',
        mode: 'markers',
        name: dataName,
        marker: { color: '#f39c12', size: 6, symbol: 'circle-open' }
    });
    
    const layout = {
        title: isUsingCustomData ? 'Custom Data & MCMC Fit' : 'Predictive Fit vs Data',
        xaxis: { title: isUsingCustomData ? 'Time' : 'Time (days)' },
        yaxis: { title: isUsingCustomData ? 'Incidence' : '% Infected' },
        height: 250,
        margin: { l: 50, r: 30, t: 30, b: 30 },
        showlegend: true,
        legend: { x: 0.7, y: 0.95, bgcolor: 'rgba(255,255,255,0.8)' }
    };
    
    Plotly.newPlot('predictivePlot', traces, layout, { responsive: true, displayModeBar: false });
}

// Plot log-likelihood trace
function plotLikelihood() {
    if (!mcmcData || mcmcData.chains.length === 0) return;
    
    const traces = [];
    const includeBurnin = document.getElementById('include-burnin')?.checked !== false;
    const burnin = includeBurnin ? 0 : parseInt(document.getElementById('burnin').value || 500);
    
    mcmcData.chains.forEach((chain, chainIndex) => {
        if (chain.logLikelihood && chain.logLikelihood.length > burnin) {
            const data = chain.logLikelihood.slice(burnin);
            traces.push({
                x: Array.from({length: data.length}, (_, i) => i + burnin),
                y: data,
                type: 'scatter',
                mode: 'lines',
                name: `Chain ${chainIndex + 1}`,
                line: { width: 1 }
            });
        }
    });
    
    const layout = {
        title: 'Log-Likelihood Trace',
        xaxis: { title: 'MCMC Step' },
        yaxis: { title: 'Log-Likelihood' },
        height: 250,
        margin: { l: 50, r: 30, t: 30, b: 30 },
        showlegend: false
    };
    
    Plotly.newPlot('likelihoodPlot', traces, layout, { responsive: true, displayModeBar: false });
}

// Calculate Rhat (potential scale reduction factor) for convergence assessment
function calculateRhat(chains, param, burnin = 0) {
    if (chains.length < 2) return NaN;
    
    const chainData = chains.map(chain => {
        if (!chain.traces[param] || chain.traces[param].length <= burnin) return [];
        return chain.traces[param].slice(burnin);
    }).filter(data => data.length > 0);
    
    if (chainData.length < 2 || chainData[0].length < 4) return NaN;
    
    const n = chainData[0].length;
    const m = chainData.length;
    
    // Calculate chain means and overall mean
    const chainMeans = chainData.map(data => data.reduce((a, b) => a + b, 0) / data.length);
    const overallMean = chainMeans.reduce((a, b) => a + b, 0) / m;
    
    // Between-chain variance
    const B = n * chainMeans.reduce((sum, mean) => sum + Math.pow(mean - overallMean, 2), 0) / (m - 1);
    
    // Within-chain variance
    const W = chainData.reduce((sum, data, i) => {
        const chainVar = data.reduce((varSum, val) => varSum + Math.pow(val - chainMeans[i], 2), 0) / (n - 1);
        return sum + chainVar;
    }, 0) / m;
    
    // Estimated variance of parameter
    const varEst = ((n - 1) / n) * W + (1 / n) * B;
    
    // Potential scale reduction factor
    const rhat = Math.sqrt(varEst / W);
    
    return rhat;
}

// Calculate effective sample size
function calculateESS(chains, param, burnin = 0) {
    if (chains.length === 0) return NaN;
    
    const chainData = chains.map(chain => {
        if (!chain.traces[param] || chain.traces[param].length <= burnin) return [];
        return chain.traces[param].slice(burnin);
    }).filter(data => data.length > 0);
    
    if (chainData.length === 0 || chainData[0].length < 4) return NaN;
    
    // Simplified ESS calculation: total samples / (1 + 2 * sum of autocorrelations)
    // For simplicity, we'll estimate autocorrelation at lag 1
    let totalSamples = 0;
    let totalAutocorr = 0;
    
    chainData.forEach(data => {
        const n = data.length;
        totalSamples += n;
        
        // Calculate lag-1 autocorrelation
        const mean = data.reduce((a, b) => a + b, 0) / n;
        let numerator = 0;
        let denominator = 0;
        
        for (let i = 0; i < n - 1; i++) {
            numerator += (data[i] - mean) * (data[i + 1] - mean);
        }
        for (let i = 0; i < n; i++) {
            denominator += Math.pow(data[i] - mean, 2);
        }
        
        const autocorr = denominator > 0 ? numerator / denominator : 0;
        totalAutocorr += Math.max(0, autocorr); // Don't let negative autocorr reduce ESS
    });
    
    const avgAutocorr = totalAutocorr / chainData.length;
    const ess = totalSamples / (1 + 2 * avgAutocorr);
    
    return Math.max(1, ess); // ESS should be at least 1
}

// Update convergence diagnostics display
function updateConvergenceDiagnostics() {
    if (!mcmcData || mcmcData.chains.length === 0) return;
    
    const includeBurnin = document.getElementById('include-burnin')?.checked !== false;
    const burnin = includeBurnin ? 0 : parseInt(document.getElementById('burnin').value || 500);
    
    let allConverged = true;
    let anyConverged = false;
    
    ['beta', 'gamma', 'sigma'].forEach(param => {
        const rhat = calculateRhat(mcmcData.chains, param, burnin);
        const ess = calculateESS(mcmcData.chains, param, burnin);
        
        const rhatStr = isNaN(rhat) ? '-' : rhat.toFixed(2);
        const essStr = isNaN(ess) ? '-' : Math.round(ess).toString();
        
        // Color coding for Rhat (good < 1.1, warning < 1.2, bad >= 1.2)
        let rhatColor = '#27ae60'; // green
        if (!isNaN(rhat)) {
            if (rhat >= 1.2) {
                rhatColor = '#e74c3c'; // red
                allConverged = false;
            } else if (rhat >= 1.1) {
                rhatColor = '#f39c12'; // orange
                allConverged = false;
            } else {
                anyConverged = true;
            }
        } else {
            allConverged = false;
        }
        
        const rhatElement = document.getElementById(`rhat-${param}`);
        const essElement = document.getElementById(`ess-${param}`);
        
        if (rhatElement) {
            rhatElement.innerHTML = `<span style="color:${rhatColor}">${rhatStr}</span>`;
        }
        if (essElement) {
            essElement.textContent = essStr;
        }
    });
    
    // Update convergence status indicator
    updateConvergenceStatus(allConverged, anyConverged);
}

// Update convergence status indicator
function updateConvergenceStatus(allConverged, anyConverged) {
    const indicatorElement = document.getElementById('convergence-indicator');
    if (!indicatorElement) return;
    
    if (allConverged) {
        indicatorElement.innerHTML = '✅ <span style="color: #27ae60;">Converged</span>';
    } else if (anyConverged) {
        indicatorElement.innerHTML = '⚠️ <span style="color: #f39c12;">Partial</span>';
    } else {
        indicatorElement.innerHTML = '❌ <span style="color: #e74c3c;">Not Converged</span>';
    }
}

// Update all MCMC plots
async function updateMCMCPlots() {
    if (mcmcData && mcmcData.currentStep > 10) {
        try {
            plotTraces();
            plotPosteriors();
            await plotPredictive();
            plotLikelihood();
            updateConvergenceDiagnostics();
        } catch (error) {
            console.error('[MCMC] Error updating plots:', error);
        }
    }
}

// Cleanup function
function cleanup() {
    if (mcmcModule) {
        // Destroy samplers and adaptive data
        for (let i = 0; i < mcmcSamplers.length; i++) {
            if (mcmcSamplers[i]) {
                mcmcModule.ccall('destroy_sampler', 'void', ['number'], [mcmcSamplers[i]]);
            }
            if (mcmcAdaptiveData[i]) {
                mcmcModule.ccall('destroy_adaptive_data', 'void', ['number'], [mcmcAdaptiveData[i]]);
            }
        }
    }
}

// Handle custom data upload
function handleDataUpload(event) {
    const file = event.target.files[0];
    if (!file) return;
    
    updateUploadStatus('Reading file...');
    
    const reader = new FileReader();
    reader.onload = function(e) {
        try {
            const text = e.target.result;
            const parsedData = parseCustomData(text);
            
            if (parsedData && parsedData.time.length > 0) {
                customData = parsedData;
                isUsingCustomData = true;
                
                // Plot the custom data
                plotCustomData(customData);
                
                // Use custom data as observed data for MCMC
                observedData = {
                    time: [...customData.time],
                    I_observed: [...customData.incidence]
                };
                
                // Enable MCMC button
                document.getElementById('mcmc-btn').disabled = false;
                updateMCMCButtonState();
                
                updateUploadStatus(`✅ Loaded ${customData.time.length} data points`);
                updateStatus(`Custom data loaded: ${customData.time.length} points, time range [${Math.min(...customData.time).toFixed(1)}, ${Math.max(...customData.time).toFixed(1)}]`);
                
                console.log('[DATA] Custom data loaded:', {
                    points: customData.time.length,
                    timeRange: [Math.min(...customData.time), Math.max(...customData.time)],
                    incidenceRange: [Math.min(...customData.incidence), Math.max(...customData.incidence)]
                });
                
            } else {
                throw new Error('No valid data found');
            }
            
        } catch (error) {
            console.error('[DATA] Upload error:', error);
            updateUploadStatus('❌ Upload failed: ' + error.message);
            customData = null;
            isUsingCustomData = false;
        }
    };
    
    reader.onerror = function() {
        updateUploadStatus('❌ Failed to read file');
        console.error('[DATA] File read error');
    };
    
    reader.readAsText(file);
}

// Parse custom data from text (CSV or space-separated)
function parseCustomData(text) {
    const lines = text.trim().split('\n');
    const time = [];
    const incidence = [];
    
    let validRows = 0;
    let errors = [];
    
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        
        // Skip empty lines and potential header lines
        if (!line || line.toLowerCase().includes('time') || line.toLowerCase().includes('incidence')) {
            continue;
        }
        
        // Try to parse as CSV (comma-separated) or space-separated
        const parts = line.includes(',') 
            ? line.split(',').map(p => p.trim()) 
            : line.split(/\s+/);
            
        if (parts.length >= 2) {
            const timeVal = parseFloat(parts[0]);
            const incidenceVal = parseFloat(parts[1]);
            
            if (!isNaN(timeVal) && !isNaN(incidenceVal) && incidenceVal >= 0) {
                time.push(timeVal);
                incidence.push(incidenceVal);
                validRows++;
            } else {
                errors.push(`Line ${i + 1}: Invalid numbers (${parts[0]}, ${parts[1]})`);
            }
        } else {
            errors.push(`Line ${i + 1}: Expected 2 columns, got ${parts.length}`);
        }
    }
    
    if (validRows < 3) {
        throw new Error(`Need at least 3 data points, found ${validRows}. Format should be [time, incidence] per line.`);
    }
    
    if (errors.length > 0 && errors.length < 5) {
        console.warn('[DATA] Parse warnings:', errors);
    }
    
    // Sort by time to ensure proper ordering
    const sortedData = time.map((t, i) => ({ time: t, incidence: incidence[i] }))
                          .sort((a, b) => a.time - b.time);
    
    return {
        time: sortedData.map(d => d.time),
        incidence: sortedData.map(d => d.incidence)
    };
}

// Plot custom data
function plotCustomData(data) {
    const trace = {
        x: data.time,
        y: data.incidence,
        type: 'scatter',
        mode: 'lines+markers',
        line: { color: '#e74c3c', width: 3 },
        marker: { color: '#e74c3c', size: 6 },
        name: 'Custom Incidence Data'
    };
    
    const layout = {
        title: 'Custom Data - Incidence over Time',
        xaxis: { title: 'Time' },
        yaxis: { title: 'Incidence' },
        showlegend: true,
        height: 400
    };
    
    Plotly.newPlot('plotContainer', [trace], layout, { responsive: true, displayModeBar: false });
}

// Update upload status display
function updateUploadStatus(message) {
    const statusElement = document.getElementById('upload-status');
    if (statusElement) {
        statusElement.textContent = message;
    }
}

// Update algorithm info display
function updateAlgorithmInfo(algorithm) {
    const descriptionElement = document.getElementById('algorithm-description');
    const referenceElement = document.getElementById('algorithm-reference');
    
    if (algorithm === 'adaptive') {
        descriptionElement.textContent = 'Auto-tuning proposal covariance';
        referenceElement.textContent = 'Andrieu & Thoms - Algorithm 4';
        referenceElement.title = 'A tutorial on adaptive MCMC';
        referenceElement.style.display = 'inline';
    } else {
        descriptionElement.textContent = 'Fixed proposal covariance';
        referenceElement.style.display = 'none';
    }
}

// Show algorithm reference modal
function showAlgorithmReference() {
    // Only show modal for adaptive algorithm
    if (mcmcAlgorithm !== 'adaptive') {
        console.log('[UI] No reference available for standard algorithm');
        return;
    }
    
    const modal = document.getElementById('referenceModal');
    const titleElement = document.getElementById('referenceTitle');
    const bodyElement = document.getElementById('referenceBody');
        titleElement.textContent = 'Adaptive Metropolis-Hastings Algorithm';
        bodyElement.innerHTML = `
            <p><strong>Reference:</strong> Andrieu, C., & Thoms, J. (2008). "A tutorial on adaptive MCMC." <em>Statistics and Computing</em>, 18(4), 343-373.</p>
            <p><strong>URL:</strong> <a href="https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf" target="_blank" style="color: #3498db;">https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf</a></p>
            <p><strong>Algorithm:</strong> Algorithm 4 (Adaptive Metropolis)</p>
            
            <h4 style="margin-top: 1.5rem; color: #2c3e50;">Algorithm Description:</h4>
            <p>The Adaptive Metropolis algorithm automatically tunes the proposal covariance matrix during sampling to improve mixing and convergence. Key features:</p>
            <ul style="margin-left: 1.5rem;">
                <li><strong>Covariance Adaptation:</strong> Uses empirical covariance from previous samples</li>
                <li><strong>Robbins-Monro Scaling:</strong> Adjusts proposal scale to target ~23.4% acceptance rate</li>
                <li><strong>Automatic Tuning:</strong> No manual tuning of proposal distributions required</li>
                <li><strong>Better Convergence:</strong> Adapts to parameter correlations and scales</li>
            </ul>
            
            <p style="margin-top: 1rem;"><strong>Implementation:</strong> This version uses a 3-phase adaptation:</p>
            <ul style="margin-left: 1.5rem;">
                <li>Steps 1-24: Fixed proposals to build sample history</li>
                <li>Steps 25-99: Robbins-Monro scaling only</li>
                <li>Steps 100+: Full covariance adaptation + scaling</li>
            </ul>
        `;
    
    modal.style.display = 'block';
}

// Hide algorithm reference modal
function hideAlgorithmReference() {
    document.getElementById('referenceModal').style.display = 'none';
}

// Close modal when clicking outside of it
window.onclick = function(event) {
    const modal = document.getElementById('referenceModal');
    if (event.target === modal) {
        hideAlgorithmReference();
    }
}

// Download posterior samples as CSV
function downloadPosteriors() {
    if (!mcmcData || !mcmcData.chains || mcmcData.chains.length === 0) {
        updateDownloadStatus('❌ No MCMC data available. Run MCMC first.');
        return;
    }
    
    const thinningFactor = parseInt(document.getElementById('thinning-factor').value);
    const includeBurnin = document.getElementById('include-burnin')?.checked !== false;
    const burnin = includeBurnin ? 0 : parseInt(document.getElementById('burnin').value || 500);
    
    // Collect all samples from all chains
    const allSamples = [];
    const chainLabels = ['Chain_1', 'Chain_2', 'Chain_3', 'Chain_4'];
    
    mcmcData.chains.forEach((chain, chainIndex) => {
        if (chain.traces.beta && chain.traces.beta.length > burnin) {
            const betaData = chain.traces.beta.slice(burnin);
            const gammaData = chain.traces.gamma.slice(burnin);
            const sigmaData = chain.traces.sigma.slice(burnin);
            
            // Apply thinning
            for (let i = 0; i < betaData.length; i += thinningFactor) {
                allSamples.push({
                    chain: chainLabels[chainIndex],
                    iteration: i + burnin,
                    beta: betaData[i],
                    gamma: gammaData[i],
                    sigma: sigmaData[i]
                });
            }
        }
    });
    
    if (allSamples.length === 0) {
        updateDownloadStatus('❌ No samples available. Check burn-in settings.');
        return;
    }
    
    // Create CSV content
    const csvHeader = 'chain,iteration,beta,gamma,sigma\n';
    const csvRows = allSamples.map(sample => 
        `${sample.chain},${sample.iteration},${sample.beta.toFixed(6)},${sample.gamma.toFixed(6)},${sample.sigma.toFixed(6)}`
    ).join('\n');
    const csvContent = csvHeader + csvRows;
    
    // Create and download file
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', `posterior_samples_thin${thinningFactor}_${new Date().toISOString().slice(0,10)}.csv`);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    
    updateDownloadStatus(`✅ Downloaded ${allSamples.length} samples (thinning: ${thinningFactor}x)`);
    console.log(`[DOWNLOAD] Exported ${allSamples.length} posterior samples with thinning factor ${thinningFactor}`);
}

// Update download status message
function updateDownloadStatus(message) {
    const statusElement = document.getElementById('download-status');
    if (statusElement) {
        statusElement.textContent = message;
    }
}

// Intervention analysis variables
let interventionData = null;
let baselineData = null;

// Waning model functions
function calculateWaningExponential(time, protection, lossRate) {
    if (lossRate === 0) {
        return protection; // No waning if loss rate is zero
    }
    return protection * Math.exp(-lossRate * time);
}

function calculateWaningGamma(time, protection, lossRate, compartments = 2) {
    if (lossRate === 0) {
        return protection; // No waning if loss rate is zero
    }
    
    // Gamma distribution with shape parameter = compartments, rate parameter = lossRate * compartments
    const shape = compartments;
    const rate = lossRate * compartments;
    
    // For gamma distribution, we need to calculate the survival function
    // P(T > t) = 1 - F(t) where F is the CDF of gamma distribution
    // Using approximation for gamma CDF
    const x = rate * time;
    let sum = 0;
    let term = 1;
    
    for (let k = 0; k < shape; k++) {
        sum += term;
        term *= x / (k + 1);
    }
    
    const cdf = 1 - Math.exp(-x) * sum;
    return protection * (1 - cdf);
}

function calculateWaningErlang3(time, protection, lossRate) {
    return calculateWaningGamma(time, protection, lossRate, 3);
}

function calculateWaningCurve(distribution, protection, lossRate, maxTime = 200) {
    const timePoints = [];
    const protectionLevels = [];
    
    for (let t = 0; t <= maxTime; t += 2) {
        timePoints.push(t);
        
        let protectionLevel;
        switch (distribution) {
            case 'exponential':
                protectionLevel = calculateWaningExponential(t, protection, lossRate);
                break;
            case 'gamma':
                protectionLevel = calculateWaningGamma(t, protection, lossRate, 2);
                break;
            case 'erlang3':
                protectionLevel = calculateWaningErlang3(t, protection, lossRate);
                break;
            default:
                protectionLevel = calculateWaningExponential(t, protection, lossRate);
        }
        
        protectionLevels.push(protectionLevel);
    }
    
    return { time: timePoints, protection: protectionLevels };
}

function plotWaningCurve() {
    const protection = parseFloat(document.getElementById('waning_protection').value) / 100;
    const lossRate = parseFloat(document.getElementById('waning_rate').value);
    const distribution = document.getElementById('waning_distribution').value;
    
    const waningData = calculateWaningCurve(distribution, protection, lossRate);
    
    const trace = {
        x: waningData.time,
        y: waningData.protection.map(p => p * 100), // Convert to percentage
        type: 'scatter',
        mode: 'lines',
        line: { color: '#e67e22', width: 2 },
        name: 'Protection Level',
        showlegend: false
    };
    
    const layout = {
        xaxis: { 
            title: '',
            range: [0, 200],
            showgrid: false,
            showticklabels: true,
            tickfont: { size: 8 },
            zeroline: false,
            automargin: true
        },
        yaxis: { 
            title: '',
            range: [0, 100],
            showgrid: false,
            showticklabels: true,
            tickfont: { size: 8 },
            zeroline: false,
            automargin: true
        },
        margin: { l: 30, r: 10, t: 10, b: 30 },
        height: 60,
        width: 140,
        showlegend: false,
        autosize: true
    };
    
    Plotly.newPlot('waningPlot', [trace], layout, { 
        responsive: true, 
        displayModeBar: false,
        staticPlot: true
    });
}

function updateWaningDistributionInfo() {
    const distribution = document.getElementById('waning_distribution').value;
    const infoElement = document.getElementById('waning_distribution_info');
    
    switch (distribution) {
        case 'exponential':
            infoElement.textContent = 'Single compartment';
            break;
        case 'gamma':
            infoElement.textContent = '2 compartments';
            break;
        case 'erlang3':
            infoElement.textContent = '3 compartments';
            break;
        default:
            infoElement.textContent = 'Single compartment';
    }
}

// Analyze intervention impact
async function analyzeIntervention() {
    if (!mcmcData || !mcmcData.bestParams || !observedData) {
        updateInterventionStatus('❌ No MCMC data available. Run MCMC first.');
        return;
    }
    
    try {
        updateInterventionStatus('🔄 Analyzing intervention impact...');
        
        // Get intervention parameters
        const coverage = parseFloat(document.getElementById('vaccine_coverage').value) / 100;
        const startTime = parseInt(document.getElementById('vaccine_start').value);
        const endTime = parseInt(document.getElementById('vaccine_end').value);
        
        // Validate parameters
        if (startTime >= endTime) {
            updateInterventionStatus('❌ Start time must be before end time');
            return;
        }
        
        if (endTime > parseInt(document.getElementById('time_period').value)) {
            updateInterventionStatus('❌ End time exceeds simulation period');
            return;
        }
        
        // Get MCMC parameters
        const beta = mcmcData.bestParams.beta;
        const gamma = mcmcData.bestParams.gamma;
        const initialInfected = parseFloat(document.getElementById('initial_infected').value);
        const timePeriod = parseInt(document.getElementById('time_period').value);
        const populationSize = parseInt(document.getElementById('population_size').value);
        
        // Run baseline simulation (no intervention)
        baselineData = await runSIRSimulation(beta, gamma, initialInfected, timePeriod, populationSize);
        
        // Run intervention simulation (with vaccination)
        interventionData = await runInterventionSimulation(
            beta, gamma, initialInfected, timePeriod, populationSize,
            coverage, startTime, endTime
        );
        
        // Calculate impact metrics
        const impactMetrics = calculateInterventionImpact(baselineData, interventionData, coverage);
        
        // Update UI
        updateInterventionResults(impactMetrics);
        plotInterventionComparison(baselineData, interventionData, startTime, endTime);
        
        // Show results section
        document.getElementById('intervention-results').style.display = 'block';
        
        // Enable download buttons
        updateInterventionDownloadButtonState();
        
        updateInterventionStatus('✅ Intervention analysis complete');
        
    } catch (error) {
        console.error('[INTERVENTION] Analysis error:', error);
        updateInterventionStatus('❌ Intervention analysis failed');
    }
}

// Run intervention simulation with waning
async function runInterventionSimulation(beta, gamma, initialInfected, timePeriod, populationSize, coverage, startTime, endTime) {
    // Get waning parameters
    const waningProtection = parseFloat(document.getElementById('waning_protection').value) / 100;
    const waningRate = parseFloat(document.getElementById('waning_rate').value);
    const waningDistribution = document.getElementById('waning_distribution').value;
    
    const time = [];
    const S = [];
    const I = [];
    const R = [];
    const V = []; // Vaccinated compartment
    
    const dt = 0.1;
    const steps = Math.floor(timePeriod / dt);
    
    // Initial conditions
    let s = 1.0 - initialInfected / 100;
    let i = initialInfected / 100;
    let r = 0.0;
    let v = 0.0; // Vaccinated
    
    for (let step = 0; step <= steps; step++) {
        const t = step * dt;
        time.push(t);
        
        // Calculate vaccination rate (linear ramp from start to end time)
        let vaccinationRate = 0;
        if (t >= startTime && t <= endTime) {
            const progress = (t - startTime) / (endTime - startTime);
            vaccinationRate = (coverage / (endTime - startTime)) * (1 - progress);
        }
        
        // Calculate current protection level based on waning model
        let currentProtection = 0;
        if (t > 0) {
            switch (waningDistribution) {
                case 'exponential':
                    currentProtection = calculateWaningExponential(t, waningProtection, waningRate);
                    break;
                case 'gamma':
                    currentProtection = calculateWaningGamma(t, waningProtection, waningRate, 2);
                    break;
                case 'erlang3':
                    currentProtection = calculateWaningErlang3(t, waningProtection, waningRate);
                    break;
                default:
                    currentProtection = calculateWaningExponential(t, waningProtection, waningRate);
            }
        }
        
        // Calculate waning rate based on distribution
        let waningToSusceptible = 0;
        
        if (waningRate === 0) {
            // No waning if loss rate is zero
            waningToSusceptible = 0;
        } else if (waningDistribution === 'exponential') {
            // Direct transition: V -> S
            waningToSusceptible = waningRate;
        } else if (waningDistribution === 'gamma') {
            // Gamma distribution with shape parameter = 2
            waningToSusceptible = waningRate * 2;
        } else if (waningDistribution === 'erlang3') {
            // Erlang-3 distribution with shape parameter = 3
            waningToSusceptible = waningRate * 3;
        }
        
        // SIRS with vaccination and waning dynamics
        // Vaccinated individuals lose immunity and become susceptible
        const dsdt = -beta * s * i - vaccinationRate * s + waningToSusceptible * v;
        const didt = beta * s * i - gamma * i;
        const drdt = gamma * i;
        
        // Vaccination dynamics: V -> S (waning)
        const dvdt = vaccinationRate * s - waningToSusceptible * v;
        
        // Update compartments
        s = Math.max(0, s + dsdt * dt);
        i = Math.max(0, i + didt * dt);
        r = Math.max(0, r + drdt * dt);
        v = Math.max(0, v + dvdt * dt);
        
        // Normalize to ensure S + I + R + V = 1
        const total = s + i + r + v;
        s /= total;
        i /= total;
        r /= total;
        v /= total;
        
        S.push(s * 100);
        I.push(i * 100);
        R.push(r * 100);
        V.push(v * 100);
    }
    
    return { time, S, I, R, V };
}

// Calculate intervention impact metrics
function calculateInterventionImpact(baseline, intervention, coverage) {
    // Find peak infections
    const baselinePeak = Math.max(...baseline.I);
    const interventionPeak = Math.max(...intervention.I);
    const peakReduction = ((baselinePeak - interventionPeak) / baselinePeak * 100).toFixed(1);
    
    // Calculate total cases averted (area under curve difference)
    let baselineCases = 0;
    let interventionCases = 0;
    
    for (let i = 0; i < baseline.I.length - 1; i++) {
        const dt = baseline.time[i + 1] - baseline.time[i];
        baselineCases += baseline.I[i] * dt;
        interventionCases += intervention.I[i] * dt;
    }
    
    const casesAverted = (baselineCases - interventionCases).toFixed(0);
    
    return {
        peakReduction: peakReduction + '%',
        casesAverted: casesAverted,
        coverageAchieved: (coverage * 100).toFixed(0) + '%'
    };
}

// Update intervention results display
function updateInterventionResults(metrics) {
    document.getElementById('peak-reduction').textContent = metrics.peakReduction;
    document.getElementById('cases-averted').textContent = metrics.casesAverted;
    document.getElementById('coverage-achieved').textContent = metrics.coverageAchieved;
}

// Plot intervention comparison with waning
function plotInterventionComparison(baseline, intervention, startTime, endTime) {
    const traces = [
        // Baseline (no intervention) - Infected only
        {
            x: baseline.time,
            y: baseline.I,
            type: 'scatter',
            mode: 'lines',
            name: 'Infected (No Intervention)',
            line: { color: '#e74c3c', width: 3, dash: 'solid' },
            yaxis: 'y'
        },
        // Intervention - Infected
        {
            x: intervention.time,
            y: intervention.I,
            type: 'scatter',
            mode: 'lines',
            name: 'Infected (With Vaccination)',
            line: { color: '#27ae60', width: 3, dash: 'solid' },
            yaxis: 'y'
        },
        // Vaccinated
        {
            x: intervention.time,
            y: intervention.V,
            type: 'scatter',
            mode: 'lines',
            name: 'Vaccinated',
            line: { color: '#9b59b6', width: 2, dash: 'solid' },
            yaxis: 'y2'
        }
    ];
    
    // Add intervention period shading
    const interventionShading = {
        x: [startTime, startTime, endTime, endTime],
        y: [0, 100, 100, 0],
        type: 'scatter',
        mode: 'lines',
        fill: 'tonexty',
        fillcolor: 'rgba(23, 162, 184, 0.1)',
        line: { color: 'rgba(23, 162, 184, 0.3)', width: 1 },
        name: 'Vaccination Period',
        showlegend: false,
        yaxis: 'y'
    };
    traces.push(interventionShading);
    
    const layout = {
        title: 'Intervention Impact: Infected vs Vaccinated',
        xaxis: { title: 'Time (days)' },
        yaxis: { 
            title: 'Infected (%)', 
            domain: [0.5, 1],
            side: 'left'
        },
        yaxis2: { 
            title: 'Vaccinated (%)', 
            domain: [0, 0.45],
            side: 'right'
        },
        height: 400,
        margin: { l: 60, r: 60, t: 40, b: 40 },
        showlegend: true,
        legend: { x: 0.7, y: 0.95, bgcolor: 'rgba(255,255,255,0.8)' }
    };
    
    Plotly.newPlot('interventionPlot', traces, layout, { responsive: true, displayModeBar: false });
}

// Plot vaccination uptake over time
function plotVaccinationUptake() {
    const coverage = parseFloat(document.getElementById('vaccine_coverage').value) / 100;
    const startTime = parseInt(document.getElementById('vaccine_start').value);
    const endTime = parseInt(document.getElementById('vaccine_end').value);
    
    // Create time points for the plot
    const timePoints = [];
    const uptakePoints = [];
    
    for (let t = 0; t <= 200; t += 2) {
        timePoints.push(t);
        
        if (t < startTime) {
            uptakePoints.push(0);
        } else if (t >= startTime && t <= endTime) {
            // Linear ramp from 0 to coverage
            const progress = (t - startTime) / (endTime - startTime);
            uptakePoints.push(coverage * progress * 100);
        } else {
            // Constant at final coverage
            uptakePoints.push(coverage * 100);
        }
    }
    
    const trace = {
        x: timePoints,
        y: uptakePoints,
        type: 'scatter',
        mode: 'lines',
        line: { color: '#17a2b8', width: 2 },
        name: 'Vaccination Uptake',
        showlegend: false
    };
    
    const layout = {
        xaxis: { 
            title: '',
            range: [0, 200],
            showgrid: false,
            showticklabels: true,
            tickfont: { size: 8 },
            zeroline: false
        },
        yaxis: { 
            title: '',
            range: [0, 100],
            showgrid: false,
            showticklabels: true,
            tickfont: { size: 8 },
            zeroline: false
        },
        margin: { l: 20, r: 10, t: 5, b: 20 },
        height: 60,
        width: 140,
        showlegend: false
    };
    
    Plotly.newPlot('uptakePlot', [trace], layout, { 
        responsive: true, 
        displayModeBar: false,
        staticPlot: true
    });
}

// Update intervention status
function updateInterventionStatus(message) {
    const statusElement = document.getElementById('intervention-status');
    if (statusElement) {
        statusElement.textContent = message;
    }
}

// Update intervention button state
function updateInterventionButtonState() {
    const interventionBtn = document.getElementById('analyze-intervention-btn');
    if (!interventionBtn) return;
    
    if (!mcmcData || !mcmcData.bestParams || !observedData) {
        interventionBtn.disabled = true;
        interventionBtn.textContent = '🔬 No MCMC Data';
        interventionBtn.style.backgroundColor = '#6c757d';
        interventionBtn.style.borderColor = '#6c757d';
    } else {
        interventionBtn.disabled = false;
        interventionBtn.textContent = '🔬 Analyze Intervention';
        interventionBtn.style.backgroundColor = '#17a2b8';
        interventionBtn.style.borderColor = '#17a2b8';
    }
    
    // Update download button states
    updateInterventionDownloadButtonState();
}

// Update intervention download button states
function updateInterventionDownloadButtonState() {
    const csvBtn = document.getElementById('download-intervention-btn');
    const plotBtn = document.getElementById('download-intervention-plot-btn');
    
    if (!csvBtn || !plotBtn) return;
    
    if (!interventionData || !baselineData) {
        csvBtn.disabled = true;
        plotBtn.disabled = true;
        csvBtn.style.backgroundColor = '#6c757d';
        csvBtn.style.borderColor = '#6c757d';
        plotBtn.style.backgroundColor = '#6c757d';
        plotBtn.style.borderColor = '#6c757d';
    } else {
        csvBtn.disabled = false;
        plotBtn.disabled = false;
        csvBtn.style.backgroundColor = '#17a2b8';
        csvBtn.style.borderColor = '#17a2b8';
        plotBtn.style.backgroundColor = '#28a745';
        plotBtn.style.borderColor = '#28a745';
    }
}

// Download intervention data as CSV
function downloadInterventionData() {
    if (!interventionData || !baselineData) {
        updateInterventionDownloadStatus('❌ No intervention data available. Run analysis first.');
        return;
    }
    
    try {
        // Prepare CSV data
        const maxLength = Math.max(baselineData.time.length, interventionData.time.length);
        const csvRows = [];
        
        // Header
        csvRows.push('time,baseline_infected,intervention_infected,intervention_vaccinated');
        
        // Data rows
        for (let i = 0; i < maxLength; i++) {
            const time = i < baselineData.time.length ? baselineData.time[i] : 
                        i < interventionData.time.length ? interventionData.time[i] : i;
            const baselineInfected = i < baselineData.I.length ? baselineData.I[i] : '';
            const interventionInfected = i < interventionData.I.length ? interventionData.I[i] : '';
            const interventionVaccinated = i < interventionData.V.length ? interventionData.V[i] : '';
            
            csvRows.push(`${time},${baselineInfected},${interventionInfected},${interventionVaccinated}`);
        }
        
        // Create and download file
        const csvContent = csvRows.join('\n');
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const link = document.createElement('a');
        const url = URL.createObjectURL(blob);
        link.setAttribute('href', url);
        link.setAttribute('download', `intervention_analysis_${new Date().toISOString().slice(0,10)}.csv`);
        link.style.visibility = 'hidden';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        
        updateInterventionDownloadStatus('✅ Intervention data downloaded successfully');
        console.log('[DOWNLOAD] Exported intervention analysis data');
        
    } catch (error) {
        console.error('[DOWNLOAD] Error exporting intervention data:', error);
        updateInterventionDownloadStatus('❌ Download failed: ' + error.message);
    }
}

// Download intervention plot as PNG
function downloadInterventionPlot() {
    if (!interventionData || !baselineData) {
        updateInterventionDownloadStatus('❌ No intervention data available. Run analysis first.');
        return;
    }
    
    try {
        // Use Plotly's download functionality
        Plotly.download('interventionPlot', 'intervention_analysis_plot.png', {
            format: 'png',
            width: 800,
            height: 600,
            scale: 2
        });
        
        updateInterventionDownloadStatus('✅ Plot downloaded successfully');
        console.log('[DOWNLOAD] Exported intervention analysis plot');
        
    } catch (error) {
        console.error('[DOWNLOAD] Error exporting plot:', error);
        updateInterventionDownloadStatus('❌ Plot download failed: ' + error.message);
    }
}

// Update intervention download status
function updateInterventionDownloadStatus(message) {
    const statusElement = document.getElementById('intervention-download-status');
    if (statusElement) {
        statusElement.textContent = message;
    }
}

// Update download button state based on MCMC data availability
function updateDownloadButtonState() {
    const downloadBtn = document.getElementById('download-posteriors-btn');
    if (!downloadBtn) return;
    
    if (!mcmcData || !mcmcData.chains || mcmcData.chains.length === 0 || 
        !mcmcData.chains[0].traces.beta || mcmcData.chains[0].traces.beta.length === 0) {
        downloadBtn.disabled = true;
        downloadBtn.textContent = '📊 No Data';
        downloadBtn.style.backgroundColor = '#6c757d';
        downloadBtn.style.borderColor = '#6c757d';
    } else {
        downloadBtn.disabled = false;
        downloadBtn.textContent = '📊 Download CSV';
        downloadBtn.style.backgroundColor = '#28a745';
        downloadBtn.style.borderColor = '#28a745';
    }
}

// Update algorithm status display in diagnostics section
function updateAlgorithmStatusDisplay(algorithm) {
    const titleElement = document.getElementById('algorithm-status-title');
    const detailElement = document.getElementById('algorithm-status-detail');
    
    if (algorithm === 'adaptive') {
        titleElement.textContent = 'Adaptive';
        titleElement.style.color = '#856404';
        detailElement.innerHTML = '1-24:Fix<br>25-99:RM<br>100+:Full';
        detailElement.style.color = '#856404';
    } else {
        titleElement.textContent = 'Standard';
        titleElement.style.color = '#6c757d';
        detailElement.innerHTML = 'Fixed<br>Covariance<br>Matrix';
        detailElement.style.color = '#6c757d';
    }
}

// Debug function to test algorithm selection
function testAlgorithmSelection() {
    const dropdown = document.getElementById('mcmc_algorithm');
    if (!dropdown) {
        console.error('[DEBUG] Algorithm dropdown not found!');
        return;
    }
    
    console.log(`[DEBUG] Current algorithm: ${mcmcAlgorithm}`);
    console.log(`[DEBUG] Dropdown value: ${dropdown.value}`);
    console.log(`[DEBUG] Available options:`, Array.from(dropdown.options).map(opt => opt.value));
    
    // Test switching
    if (mcmcAlgorithm === 'adaptive') {
        console.log('[DEBUG] Switching to standard...');
        dropdown.value = 'standard';
        dropdown.dispatchEvent(new Event('change'));
    } else {
        console.log('[DEBUG] Switching to adaptive...');  
        dropdown.value = 'adaptive';
        dropdown.dispatchEvent(new Event('change'));
    }
    
    setTimeout(() => {
        console.log(`[DEBUG] After switch - mcmcAlgorithm: ${mcmcAlgorithm}, dropdown: ${dropdown.value}`);
    }, 100);
}

// Make functions globally available for HTML onclick handlers
window.runSimulation = runSimulation;
window.startMCMC = startMCMC;
window.stopMCMC = stopMCMC;
window.handleDataUpload = handleDataUpload;
window.showAlgorithmReference = showAlgorithmReference;
window.hideAlgorithmReference = hideAlgorithmReference;
window.testAlgorithmSelection = testAlgorithmSelection;
window.downloadPosteriors = downloadPosteriors;
window.analyzeIntervention = analyzeIntervention;
window.downloadInterventionData = downloadInterventionData;
window.downloadInterventionPlot = downloadInterventionPlot;

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', initApp);

// Cleanup on page unload
window.addEventListener('beforeunload', cleanup);
