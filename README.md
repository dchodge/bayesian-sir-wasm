# Bayesian SIR WebAssembly

High-performance Bayesian SIR epidemic model with real-time adaptive MCMC parameter estimation, compiled from C++ to WebAssembly for near-native browser performance.

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
