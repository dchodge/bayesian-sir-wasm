# WebR Oscillating Sine Wave (R + C++ + WebAssembly)

A beautiful, real-time animated sine wave visualization using **WebR** (R running in WebAssembly) with a **C++ HTTP server** - completely Python-free!

## Architecture

- 🟦 **R (WebR)**: Mathematical computations and sine wave generation in the browser
- 🟩 **C++**: High-performance HTTP server with CORS support  
- 🟨 **WebAssembly**: Optional C++ math functions for enhanced performance
- 🎨 **HTML/JS**: Modern responsive UI with real-time controls

## Features

- **Pure R + C++**: No Python dependencies - truly native performance
- **Real-time Animation**: Smooth 20 FPS sine wave oscillation
- **WebR Integration**: R mathematical functions running in WebAssembly
- **C++ HTTP Server**: Lightning-fast static file serving with proper CORS
- **Interactive Controls**: Live adjustment of amplitude, frequency, speed, and phase
- **Modern UI**: Beautiful gradient design with Bootstrap styling
- **Cross-platform**: Works on macOS, Linux, and Windows

## Quick Start

### One-Command Launch
```bash
./start.sh
```
This will automatically build (if needed) and start the C++ server.

### Manual Build and Run
```bash
# Build everything
./build.sh

# Start the server
./start.sh
```

The application will automatically open in your browser at `http://localhost:8080`

## Controls

- **Amplitude**: Controls the height of the wave (0.1 - 3.0)
- **Frequency**: Controls how many cycles fit in the display (0.5 - 5.0)
- **Speed**: Controls how fast the wave moves (0.1 - 2.0)
- **Phase**: Controls the initial phase offset (0 - 6.28)

## Buttons

- **Start Animation**: Begin the real-time wave animation
- **Stop Animation**: Pause the animation
- **Reset**: Return all parameters to default values

## Build Requirements

- **CMake 3.10+**: For building C++ components
- **C++17 compiler**: GCC, Clang, or Visual Studio
- **Emscripten** (optional): For WebAssembly module compilation

### Installing Emscripten (Optional)

```bash
# Install Emscripten for WebAssembly support
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
```

## Architecture Details

### File Structure
```
├── src/
│   ├── server.cpp         # C++ HTTP Server (native)
│   ├── math.cpp           # C++ Math functions (WebAssembly)
│   └── math.hpp           # C++ Math headers
├── web/
│   ├── index.html         # Main HTML interface
│   ├── js/
│   │   └── webr-app.js   # WebR application logic
│   ├── wasm_module.js     # Generated WebAssembly wrapper
│   └── wasm_module.wasm   # Generated WebAssembly module
├── build.sh              # Build script
├── start.sh              # Launch script
├── CMakeLists.txt         # Build configuration
└── README.md             # This file
```

### Component Interaction

1. **C++ Server** serves static files with proper CORS headers
2. **WebR** loads and executes R code for sine wave mathematics
3. **WebAssembly** (optional) provides additional C++ math functions
4. **Plotly.js** renders the animated visualization
5. **JavaScript** orchestrates the real-time animation loop

### WebR Integration

The application uses WebR to execute R code like:
```r
# Generate sine wave data in R
points <- 500
x <- seq(0, 4*pi, length.out = points)
y <- amplitude * sin(frequency * x + phase + time_offset)
```

### Dependencies

- **WebR**: `https://webr.r-wasm.org/latest/webr.mjs` (CDN)
- **Plotly.js**: `https://cdn.plot.ly/plotly-latest.min.js` (CDN)
- **Bootstrap**: `https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css` (CDN)

No local dependencies required - everything loads from CDN!

## Customization

### Animation Parameters
Edit `web/js/webr-app.js`:
- `POINTS`: Number of data points (default: 500)
- `UPDATE_INTERVAL`: Animation frame delay in ms (default: 50)
- `X_RANGE`: Wave display range (default: 0 to 4π)

### Server Configuration
The C++ server accepts command-line arguments:
```bash
./build/webr_server --port 3000 --dir custom_web_dir
```

### R Code Modifications
The sine wave generation happens in R within `generateSineWaveData()`. You can modify it to create:
- Different wave types (cosine, tangent, square waves)
- Complex wave combinations
- Statistical distributions
- Custom mathematical functions

## Troubleshooting

### Build Issues
```bash
# Clean build
rm -rf build build_wasm
./build.sh
```

### WebR Not Loading
- Ensure internet connection for CDN resources
- Check browser console for errors
- Verify CORS headers (C++ server handles this automatically)

### Performance Issues
- Reduce `POINTS` in `webr-app.js`
- Increase `UPDATE_INTERVAL`
- Close other browser tabs

### Browser Compatibility
- **Required**: WebAssembly support (Chrome 57+, Firefox 52+, Safari 11+)
- **Recommended**: Modern browser with ES6 modules support

## Advanced Usage

### Custom Wave Functions
Add new wave types in the R code:
```r
# Square wave
y <- amplitude * sign(sin(frequency * x + phase + time_offset))

# Sawtooth wave  
y <- amplitude * (2 * ((frequency * x + time_offset) %% (2*pi)) / (2*pi) - 1)
```

### Multiple Waves
Combine multiple sine waves:
```r
y1 <- amplitude * sin(frequency * x + phase + time_offset)
y2 <- amplitude * 0.5 * sin(2 * frequency * x + time_offset)
y <- y1 + y2
```

## License

This project is open source and available under the MIT License.
