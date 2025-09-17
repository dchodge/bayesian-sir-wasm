// wasm-app.js - Pure C++ WebAssembly Implementation
// Simplified version with only C++ WASM support (no fallbacks)

let mcmcModule; // C++ MCMC WebAssembly module
let mcmcSamplers = []; // Array of C++ samplers for 4 chains
let mcmcAdaptiveData = []; // Array of adaptive data for 4 chains
let simulationData = null; // Store SIR simulation results
let observedData = null; // Store synthetic observed data for MCMC fitting
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
            DEFAULT_PARAMS.time_period
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
        
        // Run the SIR simulation using C++ WASM
        simulationData = await runSIRSimulation(beta, gamma, initial_infected, time_period);
        
        // Plot the results
        plotSIRResults(simulationData);
        
        // Generate synthetic observed data with noise
        observedData = generateObservedData(simulationData);
        
        // Enable MCMC button
        document.getElementById('mcmc-btn').disabled = false;
        
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
    if (mcmcRunning || !observedData) return;
    
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
    
    updateStatus('Starting MCMC sampling with C++ WebAssembly...');
    
    try {
        // Read values from sliders
        const maxSteps = parseInt(document.getElementById('mcmc_steps').value);
        const burnin = parseInt(document.getElementById('burnin').value);
        
        console.log(`[MCMC] Starting ${maxSteps} steps with ${burnin} burnin (4 adaptive chains)`);
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
                
                // Enable adaptive MCMC after adaptation start
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
                            parseInt(document.getElementById('time_period').value)
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
    }
}

// C++ Adaptive MCMC step function  
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
        
            // Call C++ adaptive MCMC step
            const success = mcmcModule.ccall('mcmc_step_adaptive', 'number', [
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
        
        if (!success) {
                throw new Error('C++ adaptive MCMC step failed');
        }
        
        // Read results
            const result = {
                params: {
            beta: mcmcModule.getValue(newBetaPtr, 'double'),
            gamma: mcmcModule.getValue(newGammaPtr, 'double'),
            sigma: mcmcModule.getValue(newSigmaPtr, 'double')
                },
                log_posterior: mcmcModule.getValue(newLogPostPtr, 'double'),
                accepted: mcmcModule.getValue(acceptedPtr, 'i32') === 1,
                used_adaptive: mcmcModule.getValue(usedAdaptivePtr, 'i32') === 1,
                acceptance_rate: mcmcModule.getValue(acceptanceRatePtr, 'double'),
                scaling_factor: mcmcModule.getValue(scalingFactorPtr, 'double')
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
    
    // Show adaptive status
    const adaptiveElement = document.getElementById('cov-status');
    if (adaptiveElement) {
        if (step > MCMC_PARAMS.adaptation_start) {
            adaptiveElement.innerHTML = `Active: 4/4 chains<br>Adaptive covariance + Robbins-Monro scaling`;
                } else {
            const stepsUntil = MCMC_PARAMS.adaptation_start - step;
            adaptiveElement.innerHTML = `Inactive<br>${stepsUntil} steps until adaptive`;
        }
    }
}

// Plot MCMC results showing original vs fitted
function plotMCMCResults(originalData, fittedData) {
    console.log('[MCMC] plotMCMCResults called with:', {
        originalDataLength: originalData?.time?.length || 0,
        fittedDataLength: fittedData?.time?.length || 0,
        originalSample: originalData?.I?.slice(0, 3),
        fittedSample: fittedData?.I?.slice(0, 3)
    });
    
    if (!originalData || !fittedData || !originalData.time || !fittedData.time) {
        console.error('[MCMC] Invalid data passed to plotMCMCResults');
            return;
        }
        
    const traces = [
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
    
    const layout = {
        title: 'MCMC Bayesian Fitting - C++ WebAssembly',
        xaxis: { title: 'Time (days)' },
        yaxis: { title: 'Percentage of Population' },
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

// Update R0 display
function updateR0Display() {
    const beta = parseFloat(document.getElementById('beta').value);
    const gamma = parseFloat(document.getElementById('gamma').value);
    const r0 = (beta / gamma).toFixed(2);
    
    document.getElementById('r0-value').textContent = r0;
    
    if (r0 > 1) {
        document.getElementById('r0-value').style.color = '#e74c3c';
        document.getElementById('epidemic-status').textContent = 'Epidemic grows';
    } else {
        document.getElementById('r0-value').style.color = '#27ae60';  
        document.getElementById('epidemic-status').textContent = 'Epidemic dies out';
    }
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
    
    // Update displays
    document.getElementById('beta-value').textContent = DEFAULT_PARAMS.beta.toFixed(2);
    document.getElementById('gamma-value').textContent = DEFAULT_PARAMS.gamma.toFixed(2);
    document.getElementById('initial_infected-value').textContent = DEFAULT_PARAMS.initial_infected.toFixed(1);
    document.getElementById('time_period-value').textContent = DEFAULT_PARAMS.time_period.toString();
    
    updateR0Display();
    
    // Run simulation with default parameters
    runSimulation();
    
    console.log('[UI] Model reset to defaults');
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
                
                // Update R0 for relevant parameters
                if (['beta', 'gamma'].includes(param.id)) {
                    updateR0Display();
                }
                
                console.log(`[UI] ${param.id} updated to ${value}`);
            });
            
            // Initialize value display
            const initialValue = parseFloat(slider.value);
            valueDisplay.textContent = initialValue.toFixed(param.decimals);
        } else {
            console.warn(`[UI] Missing slider or value display for ${param.id}`);
        }
    });
    
    // Initialize R0 display
    updateR0Display();
    
    // Add event listener for include burn-in checkbox
    const includeBurninCheckbox = document.getElementById('include-burnin');
    if (includeBurninCheckbox) {
        includeBurninCheckbox.addEventListener('change', async function() {
            console.log(`[UI] Include burn-in changed to ${this.checked}`);
            await updateMCMCPlots(); // Update plots when burn-in setting changes
        });
    }
    
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
        yaxis: { title: 'Transmission Rate (β)', domain: [0.7, 1] },
        yaxis2: { title: 'Recovery Rate (γ)', domain: [0.35, 0.65] },
        yaxis3: { title: 'Noise (σ)', domain: [0, 0.3] },
            height: 120,
            margin: { l: 45, r: 30, t: 10, b: 20 },
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
        xaxis: { title: 'Transmission Rate (β)', domain: [0, 0.32] },
        xaxis2: { title: 'Recovery Rate (γ)', domain: [0.34, 0.66] },
        xaxis3: { title: 'Noise (σ)', domain: [0.68, 1] },
        yaxis: { title: 'Density', domain: [0, 1] },
        yaxis2: { title: '', domain: [0, 1], anchor: 'x2' },
        yaxis3: { title: '', domain: [0, 1], anchor: 'x3' },
            height: 120,
            margin: { l: 45, r: 30, t: 10, b: 25 },
        showlegend: false
    };
    
    Plotly.newPlot('posteriorPlot', traces, layout, { responsive: true, displayModeBar: false });
}

// Plot predictive fit
async function plotPredictive() {
    if (!simulationData || !mcmcData.bestParams || !observedData) return;
    
    const traces = [];
    
    try {
        // Generate prediction with best parameters
        const predictedData = await runSIRSimulation(
            mcmcData.bestParams.beta,
            mcmcData.bestParams.gamma,
            parseFloat(document.getElementById('initial_infected').value),
            parseInt(document.getElementById('time_period').value)
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
    } catch (error) {
        console.warn('[MCMC] Could not generate prediction for predictive plot:', error);
    }
    
    // Add observed data points
    traces.push({
        x: observedData.time,
        y: observedData.I_observed,
        type: 'scatter',
        mode: 'markers',
        name: 'Observed',
        marker: { color: '#f39c12', size: 6, symbol: 'circle-open' }
    });
    
    const layout = {
        title: 'Predictive Fit vs Data',
        xaxis: { title: 'Time (days)' },
        yaxis: { title: '% Infected' },
            height: 120,
            margin: { l: 45, r: 30, t: 10, b: 20 },
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
            height: 120,
            margin: { l: 45, r: 30, t: 10, b: 20 },
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
    
    ['beta', 'gamma', 'sigma'].forEach(param => {
        const rhat = calculateRhat(mcmcData.chains, param, burnin);
        const ess = calculateESS(mcmcData.chains, param, burnin);
        
        const rhatStr = isNaN(rhat) ? '-' : rhat.toFixed(2);
        const essStr = isNaN(ess) ? '-' : Math.round(ess).toString();
        
        // Color coding for Rhat (good < 1.1, warning < 1.2, bad >= 1.2)
        let rhatColor = '#27ae60'; // green
        if (!isNaN(rhat)) {
            if (rhat >= 1.2) rhatColor = '#e74c3c'; // red
            else if (rhat >= 1.1) rhatColor = '#f39c12'; // orange
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

// Make functions globally available for HTML onclick handlers
window.runSimulation = runSimulation;
window.startMCMC = startMCMC;
window.stopMCMC = stopMCMC;

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', initApp);

// Cleanup on page unload
window.addEventListener('beforeunload', cleanup);
