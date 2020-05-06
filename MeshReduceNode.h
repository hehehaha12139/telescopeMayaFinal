#pragma once
#include <maya/MPxNode.h>
#include <maya/MStatus.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MDataBlock.h>
#include <maya/MObject.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include "mesh.h"
#include <queue>
#include "curve.h"
#include <unordered_set>

#define McheckErr(stat,msg)                     \
        if ( MS::kSuccess != stat ) {   \
                cerr << msg;                            \
                return MS::kFailure;            \
        }

struct BulbPoints 
{
	bool isFirst;
	int index;
	float radius;
	glm::vec3 pos;
	int count = 0;
};

struct IntPair 
{
	int first;
	int second;

	IntPair(int first, int second) : first(first), second(second) {};
};

class MeshReduceNode : public MPxNode
{
public:
	MeshReduceNode() {};
	virtual ~MeshReduceNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static void* creator();
	static MStatus initialize();
	static bool hasComputed;
	std::vector<glm::mat4> vertQuadric;
	std::vector<glm::mat4> faceQuadric;
	std::vector<float> edgeCosts;
	std::priority_queue<FaceNode*, std::vector<FaceNode*>, cmp> PQ;
	std::vector<std::unique_ptr<FaceNode>> faceNodes;
	std::vector<std::unordered_set<int>*> adjacency;
	std::vector<std::unique_ptr<BulbPoints>> outPoints;
	std::vector<std::vector<glm::vec3>> outSplineKeys;
	std::vector<std::unique_ptr<std::unordered_set<int>>> adjStore;
	std::vector<std::unique_ptr<Curve>> discreteCurves;
	std::vector<std::unique_ptr<Cube>> otherCubes;
	std::vector<std::vector<unique_ptr<Cube>>> cubes;
	UnionFind uf;

	static MTypeId id;

	Mesh halfEdgeMesh;

	// Telescope Attributes
	//static MObject impulseCount;
	static MObject inputFileName;
	static MObject outputCurves;
	static MObject outputJoints;
	static MObject outputCurveIndices;
	static MObject outputCurveCount;
	static MObject outputJointCount;
	static MObject outputCtrlPointCount;
	static MObject outputGeo;
	static MObject outputTelGeo;
	static bool enterOnce;

	float edgeQuadric(HalfEdge* he);
	
	void FixPQ(Vertex* v);
	void InitQuadrics();
	void InitPQ();
	void CollapseAll();

	HalfEdge* optimalEdge(Face* f);
	HalfEdge* optimalEdgeShortest(Face* f);
	float Cost(HalfEdge* he);
	bool isBoundary(Face* f);
	float lengthCost(HalfEdge* he);
};
