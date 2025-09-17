#!/bin/bash

# Setup Emscripten for C++ WebAssembly compilation
echo "🔧 Setting up Emscripten for C++ WebAssembly..."

# Check if emsdk already exists
if [ -d "emsdk" ]; then
    echo "✅ Emscripten SDK already exists"
    cd emsdk
    source ./emsdk_env.sh
    cd ..
else
    echo "📦 Downloading Emscripten SDK..."
    git clone https://github.com/emscripten-core/emsdk.git
    cd emsdk
    
    echo "🔧 Installing latest Emscripten..."
    ./emsdk install latest
    ./emsdk activate latest
    
    echo "🌐 Setting up environment..."
    source ./emsdk_env.sh
    cd ..
    
    echo "✅ Emscripten installed successfully!"
fi

echo ""
echo "🚀 Building C++ WebAssembly module..."

# Create build directory for WASM
mkdir -p build_wasm
cd build_wasm

# Build with Emscripten
emcmake cmake .. -DCMAKE_BUILD_TYPE=Release
if [ $? -ne 0 ]; then
    echo "❌ CMake configuration failed"
    exit 1
fi

make -j$(nproc 2>/dev/null || echo 4)
if [ $? -ne 0 ]; then
    echo "❌ Build failed"
    exit 1
fi

cd ..

# Copy WASM files to web directory
echo "📁 Copying WebAssembly files..."
cp build_wasm/web/wasm_module.js web/ 2>/dev/null
cp build_wasm/web/wasm_module.wasm web/ 2>/dev/null

if [ -f "web/wasm_module.js" ] && [ -f "web/wasm_module.wasm" ]; then
    echo "✅ C++ WebAssembly module built successfully!"
    echo "📂 Files: web/wasm_module.js, web/wasm_module.wasm"
    echo ""
    echo "🌊 Your sine wave now has 3-tier fallback:"
    echo "   1. 🟦 WebR (R in browser)"
    echo "   2. 🟨 C++ WebAssembly (your code!)"
    echo "   3. 🟩 JavaScript (final fallback)"
else
    echo "⚠️  WebAssembly files not found, but build completed"
fi

echo ""
echo "🚀 Ready! Start your server with: ./start.sh"
