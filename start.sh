#!/bin/bash

# WebR Oscillating Sine Wave - Quick Start Script (R + C++ + WebAssembly)
echo "🌊 Starting WebR Oscillating Sine Wave Server..."
echo ""

# Check if C++ server is built
if [ ! -f "build/webr_server" ]; then
    echo "🔧 C++ server not found. Building first..."
    
    # Try CMake build first, fall back to simple build
    if command -v cmake &> /dev/null; then
        ./build.sh
        BUILD_RESULT=$?
    else
        echo "⚠️  CMake not found, using direct compilation..."
        ./build_simple.sh
        BUILD_RESULT=$?
    fi
    
    if [ $BUILD_RESULT -ne 0 ]; then
        echo "❌ Build failed. Please fix any errors and try again."
        exit 1
    fi
fi

# Check if web directory exists
if [ ! -d "web" ]; then
    echo "❌ Web directory not found!"
    echo "Please make sure the web directory exists with your HTML files."
    exit 1
fi

echo "🚀 Starting C++ HTTP Server..."
echo "📁 Web root: web/"
echo "🔗 URL: http://localhost:8080"
echo "⚡ Features: R (webR) + C++ + WebAssembly"
echo ""

# Start the C++ server
./build/webr_server --port 8080 --dir web

