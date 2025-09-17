#!/bin/bash

# Simple Direct C++ Compilation (no CMake required)
echo "🌊 Building WebR Oscillating Sine Wave Application (Direct Compilation)"
echo ""

# Check if we have a C++ compiler
if command -v g++ &> /dev/null; then
    CXX=g++
elif command -v clang++ &> /dev/null; then
    CXX=clang++
elif command -v c++ &> /dev/null; then
    CXX=c++
else
    echo "❌ No C++ compiler found. Please install one of:"
    echo "   - GCC (g++)"
    echo "   - Clang (clang++)"
    echo "   - Or install CMake to use the full build system"
    exit 1
fi

echo "✅ Using C++ compiler: $CXX"

# Create build directory
mkdir -p build

# Build C++ HTTP Server directly
echo "🚀 Compiling C++ HTTP Server..."

$CXX -std=c++17 -O2 -pthread \
    src/server.cpp \
    -o build/webr_server

if [ $? -ne 0 ]; then
    echo "❌ Failed to compile C++ server"
    exit 1
fi

echo "✅ C++ Server compiled successfully!"
echo ""

# Note about WebAssembly
echo "📝 Note: WebAssembly module requires Emscripten"
echo "   The WebR visualization will work without it"
echo "   To install Emscripten: https://emscripten.org"
echo ""

echo "🎉 Build complete!"
echo "✅ C++ HTTP Server: build/webr_server"
echo ""
echo "🚀 To start the server, run: ./start.sh"
