#!/bin/bash

# WebR Oscillating Sine Wave - Build Script (R + C++ + WebAssembly)
echo "🌊 Building WebR Oscillating Sine Wave Application"
echo ""

# Clean build directory
if [ -d "build" ]; then
    echo "🧹 Cleaning previous build..."
    rm -rf build
fi

mkdir -p build
cd build

# Build C++ HTTP Server (native)
echo "🚀 Building C++ HTTP Server..."
cmake .. -DCMAKE_BUILD_TYPE=Release
if [ $? -ne 0 ]; then
    echo "❌ Failed to configure C++ server build"
    exit 1
fi

make -j$(nproc)
if [ $? -ne 0 ]; then
    echo "❌ Failed to build C++ server"
    exit 1
fi

echo "✅ C++ Server built successfully"

# Check if Emscripten is available for WebAssembly build
if command -v emcc &> /dev/null; then
    echo "🌊 Building WebAssembly module..."
    
    # Clean and rebuild with Emscripten
    cd ..
    if [ -d "build_wasm" ]; then
        rm -rf build_wasm
    fi
    mkdir -p build_wasm
    cd build_wasm
    
    # Build WASM module
    emcmake cmake .. -DCMAKE_BUILD_TYPE=Release
    if [ $? -ne 0 ]; then
        echo "❌ Failed to configure WebAssembly build"
        exit 1
    fi
    
    make -j$(nproc)
    if [ $? -ne 0 ]; then
        echo "❌ Failed to build WebAssembly module"
        exit 1
    fi
    
    echo "✅ WebAssembly module built successfully"
    
    # Copy WASM files to web directory
    cd ..
    cp build_wasm/web/* web/ 2>/dev/null || true
    
else
    echo "⚠️  Emscripten not found - skipping WebAssembly build"
    echo "   The WebR visualization will still work without the C++ WASM module"
    echo "   To install Emscripten: https://emscripten.org/docs/getting_started/downloads.html"
fi

echo ""
echo "🎉 Build complete!"
echo "✅ C++ HTTP Server: build/webr_server"
if command -v emcc &> /dev/null; then
    echo "✅ WebAssembly Module: web/wasm_module.js, web/wasm_module.wasm"
fi
echo ""
echo "🚀 To start the server, run: ./start.sh"
