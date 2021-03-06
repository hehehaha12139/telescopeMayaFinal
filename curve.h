#ifndef CURVE_H
#define CURVE_H

#include "discrete.h"
#include "cubearray.h"
#include "cube.h"
#include <maya/MPointArray.h>
#include <maya/MPoint.h>
#include <maya/MIntArray.h>
class Shell
{
public:
    Shell();

    void addCylinder(const std::vector<std::vector<glm::vec4>>& cylinder, glm::vec4 color = _red);

    vector<int> ib;
    vector<glm::vec4> vb;

    vector<vector<vector<glm::vec4>>> cylinders;

    glm::mat4 transform, animatedTransform, junctureAnimatedTransform;
    MPointArray points;
    MIntArray faceConnects;
    MIntArray faceCounts;
};

class Curve
{
public:
    Curve(Cube* parentCube, Cube* childCube);
    Curve();
    vector<glm::vec4>* points;
    vector<DCurvePoint> discretePoints;
    vector<glm::vec4> torsionImpulsePoints;
    glm::vec3 startingPoint;
    glm::vec3 startingTangent;
    glm::vec3 startingBinormal;
    glm::vec3 targetEndPoint;
    bool computed;
    bool processed;

    // Compute discrete points
    void discretilize(float segLength = 0.1f);
    // Calculte shell fragment number of given curve
    int computeNumImpulses();

    static int numImpulses;

    vector<float> evenlySpacePoints(int num);

    float arcLength;

    void makeImpulseCurve();
    void makeTelescope(float radius = 0.8f);
    void makeShells();

    vector<vector<glm::vec4>> generateCylinder(TelescopeParameters tParams, float nextRadius = 0.f);
    vector<glm::vec4> generateCircle(int circNum, glm::vec3 centerPoint, glm::vec3 direction,
        glm::vec3 normal, float radius);
    int VERTS_PER_CIRCLE = 30;//100;
    int CUTS_PER_CYLINDER = 10;//40;
    glm::fquat getLocalRotationAlongPath(float t, float curvature, float torsion, float length);
    vector<unique_ptr<Shell>> shells;

    glm::mat4 getCurrentTransform(TelescopeParameters tParams, float t);

    enum ExtensionState { RETRACTED, EXTENDED, RETRACTING, EXTENDING };
    ExtensionState extensionState;
    float extensionExtent;

    vector<TelescopeParameters> tParams;

    float segmentLength;

    void reAssignPoints();

    unique_ptr<vector<CurveSegment>> pSegments;

    float calcArcLength();

    void AddPointsOfSegment(CurveSegment seg);
    glm::vec3 transformedHelixPoint(CurveSegment cs, float arcLen);

    bool hasAssigned, hasTelescope;


    glm::vec3 childBasedPosition(CurveSegment parent, CurveSegment child);
    glm::mat3 childBasedRotation(CurveSegment parent, CurveSegment child);
    void updateSegmentTransforms();



    std::vector<glm::vec4> keys, ctrlPoints, curve;
    void updateCurve();
    // interpolation
    void generateKeys();
    void computeCtrlPoints();
    void interpolate();
    glm::vec4 interpolateSegment(int segment, float t);
    Cube* parentCube, * childCube;
    vector<Cube*> curveCubes;
    vector<Curve*> childrenCurves;
    float endRadius;
};



#endif // CURVE_H

