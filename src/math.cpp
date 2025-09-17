// src/math.cpp
#include "math.hpp"
#include <cmath>

std::vector<double> FigureGenerator::generateSineWave(int points, double amplitude, double frequency) {
    std::vector<double> result;
    for (int i = 0; i < points; i++) {
        double x = (2 * M_PI * i) / points;
        double y = amplitude * sin(frequency * x);
        result.push_back(y);
    }
    return result;
}

using namespace emscripten;

EMSCRIPTEN_BINDINGS(module) {
    class_<FigureGenerator>("FigureGenerator")
        .constructor<>()
        .function("generateSineWave", &FigureGenerator::generateSineWave);
    register_vector<double>("vector<double>");
}