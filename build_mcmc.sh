#!/bin/bash

echo "🔥 Building C++ MCMC WebAssembly Module..."

# Check if emscripten is available
if ! command -v emcc &> /dev/null; then
    echo "📦 Setting up Emscripten..."
    
    # Check if emsdk directory exists
    if [ ! -d "emsdk" ]; then
        echo "⬇️  Downloading Emscripten SDK..."
        git clone https://github.com/emscripten-core/emsdk.git
    fi
    
    cd emsdk
    echo "🔧 Installing latest Emscripten..."
    ./emsdk install latest
    ./emsdk activate latest
    source ./emsdk_env.sh
    cd ..
fi

echo "🚀 Compiling C++ MCMC to WebAssembly..."

# Compile with optimizations and exports
emcc src/mcmc.cpp -o web/mcmc_module.js \
  -s EXPORTED_FUNCTIONS="['_create_sampler', '_destroy_sampler', '_create_adaptive_data', '_destroy_adaptive_data', '_reset_adaptive_data', '_run_sir_simulation', '_mcmc_step', '_mcmc_step_adaptive', '_enable_adaptive_sampling', '_update_empirical_stats', '_add_sample_to_adaptive', '_get_adaptive_info', '_malloc', '_free']" \
  -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap', 'getValue', 'setValue']" \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s MODULARIZE=1 \
  -s EXPORT_NAME="'MCMCModule'" \
  -O3 \
  -std=c++17 \
  --bind

if [ $? -eq 0 ]; then
    echo "✅ C++ MCMC WebAssembly module built successfully!"
    echo "📁 Output: web/mcmc_module.js + web/mcmc_module.wasm"
    ls -la web/mcmc_module.*
else
    echo "❌ Build failed!"
    exit 1
fi
