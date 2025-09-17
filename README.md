# Bayesian SIR WebAssembly

High-performance Bayesian SIR epidemic model with real-time adaptive MCMC parameter estimation, compiled from C++ to WebAssembly for near-native browser performance.

[![Build and Deploy](https://github.com/YOUR_USERNAME/sir_bayes/actions/workflows/deploy.yml/badge.svg)](https://github.com/YOUR_USERNAME/sir_bayes/actions/workflows/deploy.yml)

## 🚀 Live Demo

**[Try the application here →](https://YOUR_USERNAME.github.io/sir_bayes/)**

## Overview

A fast, interactive web application for Bayesian inference of SIR (Susceptible-Infected-Recovered) epidemic models using advanced adaptive MCMC algorithms. Built with C++ and compiled to WebAssembly for 10-50x performance improvements over pure JavaScript implementations.

## Key Features

- ⚡ **High Performance**: C++ WebAssembly core with near-native speed
- 🧠 **Advanced MCMC**: Adaptive algorithms with automatic parameter tuning
- 📊 **Real-time Visualization**: Interactive parameter traces and posterior distributions
- 📈 **Custom Data Import**: Upload your own [time, incidence] data for analysis  
- 🎯 **Bayesian Inference**: Full uncertainty quantification with convergence diagnostics
- 🔬 **Research-Ready**: Publication-quality plots and comprehensive diagnostics

## Quick Start

### One-Command Launch
```bash
./start.sh
```

### Manual Build
```bash
# Build WebAssembly module
./build_mcmc.sh

# Build C++ server
./build_simple.sh

# Start server
./start.sh
```

## Using Custom Data

The application supports loading your own epidemic data for analysis:

1. **Format**: Prepare data as CSV or TXT with two columns: `[time, incidence]`
2. **Upload**: In Section 1, click "Upload [time, incidence]" to load your data
3. **Analyze**: The app will plot your data and enable MCMC fitting
4. **Sample**: Download `sample_data.csv` for format reference

### Data Format Example
```csv
time,incidence
0.0,1.2
5.0,2.1
10.0,4.8
15.0,8.9
```

The system automatically validates data, handles both comma and space-separated formats, and sorts by time for consistent analysis.

## MCMC Algorithms

### Adaptive Metropolis-Hastings
- **Reference**: [Andrieu & Thoms (2008)](https://people.eecs.berkeley.edu/~jordan/sail/readings/andrieu-thoms.pdf) - Algorithm 4
- **Features**: Automatic covariance adaptation and Robbins-Monro scaling
- **Benefits**: No manual tuning, better convergence, adapts to parameter correlations

### Standard Metropolis-Hastings  
- **Reference**: Metropolis et al. (1953)
- **Features**: Fixed proposal distributions with pre-tuned parameters
- **Benefits**: Simple, well-understood, reliable baseline performance


This project is licensed under the MIT License - see the LICENSE file for details.


## Acknowledgments

- **MCMC Theory**: Andrieu & Thoms (2008), Haario et al. (2001)
- **WebAssembly**: Emscripten project
- **Visualization**: Plotly.js team
- **Mathematical Foundation**: Metropolis et al. (1953)
