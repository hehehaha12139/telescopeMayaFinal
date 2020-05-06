#pragma once
#include "cube.h"

Cube::Cube()
    : transform(glm::mat4(0.0f)), type(NORMAL), selected(false), rootCube(nullptr), parentJuncture(nullptr), parentCurve(nullptr), neighbor(nullptr), junctureId(-1),processed(false)
{}

Cube::Cube(glm::mat4 transform, int type)
    : transform(transform), selected(false), type(type), rootCube(nullptr), parentJuncture(nullptr), parentCurve(nullptr), neighbor(nullptr), junctureId(-1),processed(false)
{}

