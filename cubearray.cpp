#include "cubearray.h"

void CubeArray::computeCtrlPoints()
{
    ctrlPoints.clear();
    if (keys.size() <= 1) return;
    glm::vec4 startPoint = keys[0] + 0.25f * (keys[0] - keys[1]);
    glm::vec4 endPoint = keys[keys.size() - 1] + 0.25f * (keys[keys.size() - 1] - keys[keys.size() - 2]);
    for (int i = 1; i < keys.size(); ++i)
    {
        glm::vec4 b0, b1, b2, b3, t0, t3;
        b0 = keys[i - 1], b3 = keys[i];
        if (i == 1) t0 = (keys[i] - startPoint) / 2.f;
        else t0 = (keys[i] - keys[i - 2]) / 2.f;
        b1 = b0 + t0 / 3.f;
        if (i == keys.size() - 1) t3 = (endPoint - keys[i - 1]) / 2.f;
        else t3 = (keys[i + 1] - keys[i - 1]) / 2.f;
        b2 = b3 - t3 / 3.f;
        std::vector<glm::vec4> v{b0, b1, b2, b3};
        ctrlPoints.insert(ctrlPoints.end(), v.begin(), v.end());
    }
}

void CubeArray::interpolate()
{
    curve.clear();
    if (keys.size() <= 1) return;
    for (int segment = 0; segment < keys.size() - 1; ++segment)
    {
        for (float t = 0.f; t < 1.f - std::numeric_limits<float>::min(); t += 0.01f)
        {
            curve.push_back(interpolateSegment(segment, t));
        }
    }
    curve.push_back(interpolateSegment(keys.size() - 2, 1.f));
}

glm::vec4 CubeArray::interpolateSegment(int segment, float t)
{
    glm::vec4 b0 = ctrlPoints[segment * 4 + 0],
            b1 = ctrlPoints[segment * 4 + 1],
            b2 = ctrlPoints[segment * 4 + 2],
            b3 = ctrlPoints[segment * 4 + 3],
            p00 = b0 + (b1 - b0) * t,
            p10 = b1 + (b2 - b1) * t,
            p20 = b2 + (b3 - b2) * t,
            p01 = p00 + (p10 - p00) * t,
            p11 = p10 + (p20 - p10) * t,
            p02 = p01 + (p11 - p01) * t;
    return p02;
}

void CubeArray::updateCurve()
{
    computeCtrlPoints();
    interpolate();
}
