# CMakeLists.txt - WebR Oscillating Sine Wave (R + C++ + WebAssembly)
cmake_minimum_required(VERSION 3.10)
project(WebRSineWave)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source directory
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Check if we're building with Emscripten (for WASM) or native (for server)
if(EMSCRIPTEN)
    message(STATUS "🌊 Building WebAssembly module with Emscripten")
    
    # WebAssembly module
    add_executable(wasm_module "${SOURCE_DIR}/math.cpp")
    
    # Enable embind for this target
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --bind")
    
    # Set output directory
    set_target_properties(wasm_module PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/web"
    )
    
    # Set Emscripten link flags
    set_target_properties(wasm_module PROPERTIES
        LINK_FLAGS "-s WASM=1 \
                    -s 'EXPORTED_RUNTIME_METHODS=[\"ccall\",\"cwrap\"]' \
                    -s 'EXPORTED_FUNCTIONS=[\"_malloc\",\"_free\"]' \
                    -s ALLOW_MEMORY_GROWTH=1 \
                    -s NO_EXIT_RUNTIME=1 \
                    -s MODULARIZE=1 \
                    -s 'EXPORT_NAME=\"createModule\"' \
                    --no-entry"
    )
    
else()
    message(STATUS "🚀 Building native C++ HTTP server")
    
    # Find required packages for the server
    find_package(Threads REQUIRED)
    
    # C++ HTTP Server
    add_executable(webr_server "${SOURCE_DIR}/server.cpp")
    
    # Link threads library
    target_link_libraries(webr_server PRIVATE Threads::Threads)
    
    # Set output directory
    set_target_properties(webr_server PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
    
    # Enable filesystem library support
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
        target_link_libraries(webr_server PRIVATE stdc++fs)
    endif()
    
endif()

# Custom target to build both WASM and server
add_custom_target(build_all 
    COMMENT "Building both WebAssembly module and C++ server"
)