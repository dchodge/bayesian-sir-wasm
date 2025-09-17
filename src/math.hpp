#pragma once
#include <vector>
#include <emscripten/bind.h>

class FigureGenerator {
public:
    FigureGenerator() = default;
    std::vector<double> generateSineWave(int points, double amplitude, double frequency);
};