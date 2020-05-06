#pragma once
#ifndef CUBE_H
#define CUBE_H
#define SIZE 0.5f
#include <glm/glm.hpp>

class Curve;
class Cube
{
public:
    Cube();
    Cube(glm::mat4 transform, int type = NORMAL);

    glm::mat4 transform;

    bool selected;
    int type;

    Cube* rootCube, * parentJuncture;
    Curve* parentCurve;
    Cube* neighbor;
    int junctureId;
    bool processed;
    enum CUBETYPE { RED, GREEN, BLUE, NORMAL, GENERATOR, JUNCTURE };

};
#endif // CUBE_H
