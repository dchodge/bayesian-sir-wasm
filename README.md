# An Interactive Bayesian Epidemic Modeling Platform

A web-based platform for Bayesian inference of compartmental epidemic models with real-time parameter estimation, vaccine intervention analysis, and waning immunity modeling.


[![Build and Deploy](https://github.com/dchodge/bayesian-sir-wasm/actions/workflows/deploy.yml/badge.svg)](https://github.com/dchodge/bayesian-sir-wasm/actions/workflows/deploy.yml)

## 🚀 Live Demo

**[Try the application here →](https://dchodge.github.io/bayesian-sir-wasm/)**

## What is this?

It is an interactive web application that enables users to perform Bayesian analysis of epidemic data using compartmental models. The platform combines statistical methods with high-performance computing to provide real-time parameter estimation, uncertainty quantification, and intervention analysis.

### Mathematical Foundation

**Compartmental Models**: The platform implements SIR (Susceptible-Infected-Recovered) and extended SIRS models with vaccination compartments. The core differential equations are:

```
dS/dt = -βSI - vS + γR + wV
dI/dt = βSI - γI  
dR/dt = γI - γR
dV/dt = vS - wV
```

Where:
- **β**: Transmission rate (infectiousness)
- **γ**: Recovery rate (1/infectious period)
- **v**: Vaccination rate
- **w**: Waning immunity rate
- **S, I, R, V**: Population compartments (susceptible, infected, recovered, vaccinated)

**Vaccine Waning Models**: The platform supports three waning distributions:
- **Exponential**: Simple exponential decay (1 compartment)
- **Gamma**: Gamma distribution with shape parameter = 2 (2 compartments)  
- **Erlang-3**: Erlang-3 distribution with shape parameter = 3 (3 compartments)

### Statistical Methods

**Bayesian Inference**: The platform uses Markov Chain Monte Carlo (MCMC) methods to estimate model parameters and quantify uncertainty:

**Likelihood Function**:
```
L(β,γ,σ) = ∏ᵢ N(yᵢ|I(tᵢ),σ²)
```

Where yᵢ are observed incidence data and I(tᵢ) are model predictions.

**MCMC Algorithms**:
- **Adaptive Metropolis-Hastings**: Automatic covariance adaptation and Robbins-Monro scaling
- **Standard Metropolis-Hastings**: Fixed proposal distributions for baseline comparison
- **Convergence Diagnostics**: Gelman-Rubin R̂ statistic and Effective Sample Size (ESS)

**Prior Distributions**: Fixed prior specification for all parameters with automatic scaling.

### Software Architecture

**High-Performance Computing**: 
- **C++ Core**: Mathematical computations and MCMC algorithms
- **WebAssembly Compilation**: Near-native browser performance (10-50x faster than JavaScript)
- **Emscripten Toolchain**: C++ to WASM compilation

**Web Technologies**:
- **Frontend**: HTML5, CSS3, JavaScript (ES6+)
- **Visualization**: Plotly.js for interactive charts and real-time updates
- **Data Handling**: CSV/TXT import/export with automatic validation

## Key Features

- ⚡ **High Performance**: C++ WebAssembly core with near-native speed (10-50x faster than JavaScript)
- 🧠 **Advanced MCMC**: Adaptive Metropolis-Hastings with automatic parameter tuning and convergence diagnostics
- 📊 **Real-time Visualization**: Interactive parameter traces, posterior distributions, and epidemic curves
- 📈 **Custom Data Import**: Upload your own [time, incidence] data for analysis with automatic validation
- 🎯 **Bayesian Inference**: Full uncertainty quantification with Gelman-Rubin R̂ and Effective Sample Size
- 💉 **Vaccine Modeling**: Comprehensive vaccination intervention analysis with waning immunity
- 📉 **Waning Distributions**: Exponential, Gamma, and Erlang-3 models for vaccine waning
- 🔬 **Research-Ready**: Publication-quality plots, comprehensive diagnostics, and data export
- 📱 **Interactive Interface**: Real-time parameter adjustment with immediate visual feedback
- 📊 **Data Export**: Download simulation results as CSV and plots as PNG for further analysis

## How It Works

### User Workflow

1. **Data Input**: Upload your epidemic data (time, incidence) or use the built-in sample data
2. **Parameter Setup**: Adjust initial conditions, population size, and model parameters
3. **MCMC Analysis**: Run Bayesian inference to estimate transmission and recovery rates
4. **Convergence Monitoring**: Real-time diagnostics show when chains have converged
5. **Intervention Analysis**: Model vaccination scenarios with different waning assumptions
6. **Results Export**: Download data and plots for further analysis or publication

### Technical Workflow

1. **Data Validation**: Automatic format checking and time-series sorting
2. **Model Compilation**: C++ code compiled to WebAssembly for browser execution
3. **MCMC Sampling**: Adaptive algorithms run in parallel with real-time updates
4. **Convergence Assessment**: Gelman-Rubin R̂ and ESS calculated continuously
5. **Visualization**: Plotly.js renders interactive charts with live updates
6. **Export Generation**: CSV and PNG files created client-side for download

### Performance Characteristics

- **MCMC Speed**: 1000+ iterations per second on modern hardware
- **Memory Usage**: Efficient C++ memory management with automatic garbage collection
- **Browser Compatibility**: Works in all modern browsers supporting WebAssembly
- **Data Limits**: Handles datasets with thousands of time points efficiently

## Acknowledgments

- **MCMC Theory**: Andrieu & Thoms (2008), Haario et al. (2001)
- **WebAssembly**: Emscripten project
- **Visualization**: Plotly.js team
- **Mathematical Foundation**: Metropolis et al. (1953)
