#ifndef CUBEARRAY_H
#define CUBEARRAY_H

#include <vector>
#include <glm/glm.hpp>

class CubeArray
{
public:
    std::vector<glm::vec4> keys, ctrlPoints, curve;
    void updateCurve();
    void computeCtrlPoints();
    void interpolate();
    glm::vec4 interpolateSegment(int segment, float t);
};

#endif // CUBEARRAY_H
