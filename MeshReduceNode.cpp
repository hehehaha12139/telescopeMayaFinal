#include "MeshReduceNode.h"
#include <glm/glm.hpp>
#include <maya/MTime.h>
#include <maya/MFnMeshData.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MFnArrayAttrsData.h>
#include <maya/MVectorArray.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MItMeshVertex.h>
#include "cubeArray.h"
#include "global.h"
#include <queue>
#include <stack>
#include <map>


MTypeId MeshReduceNode::id(0x80003);
MObject MeshReduceNode::inputFileName;
MObject MeshReduceNode::outputCurves;
MObject MeshReduceNode::outputJoints;
MObject MeshReduceNode::outputCurveIndices;
MObject MeshReduceNode::outputCurveCount;
MObject MeshReduceNode::outputJointCount;
MObject MeshReduceNode::outputCtrlPointCount;
MObject MeshReduceNode::outputGeo;
MObject MeshReduceNode::outputTelGeo;

bool MeshReduceNode::hasComputed = false;
bool MeshReduceNode::enterOnce = true;

void* MeshReduceNode::creator()
{
	return new MeshReduceNode;
}

MStatus MeshReduceNode::initialize()
{
	MFnNumericAttribute numAttr;
	MFnUnitAttribute unitAttr;
	MFnTypedAttribute typeAttr;
	hasComputed = false;
	MStatus returnStatus;

	cout << "InitializeNode!" << endl;

	MVectorArray defaultPoints;
	MFnVectorArrayData defaultArray;
	MObject defaultAttr;

	defaultPoints.clear();
	defaultAttr = defaultArray.create(defaultPoints);

	MeshReduceNode::inputFileName = typeAttr.create("inFileName", "ifn", MFnData::kString, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve keys attribute");

	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	MeshReduceNode::outputCurveCount = numAttr.create("outCurveCount", "occ", MFnNumericData::kInt, 0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setStorable(false);

	MeshReduceNode::outputJointCount = numAttr.create("outJointCount", "ojc", MFnNumericData::kInt, 0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setStorable(false);

	MeshReduceNode::outputCtrlPointCount = numAttr.create("outCtrlPointCount", "ocp", MFnNumericData::kInt, 0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setStorable(false);


	MeshReduceNode::outputCurves = numAttr.create("outCurves", "oc", MFnNumericData::k4Double, 0.0f, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setArray(true);
	numAttr.setStorable(false);

	MeshReduceNode::outputJoints = numAttr.create("outJoints", "oj", MFnNumericData::k3Double, 0.0f, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setArray(true);
	numAttr.setStorable(false);

	MeshReduceNode::outputCurveIndices = typeAttr.create("outIndices", "od", MFnData::kIntArray, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setStorable(false);

	MeshReduceNode::outputGeo = typeAttr.create("outGeo", "oGe", MFnData::kMesh, &returnStatus);
	McheckErr(returnStatus, "ERROR creating lSystem output geometry attribute\n");
	numAttr.setStorable(false);

	MeshReduceNode::outputTelGeo = typeAttr.create("outTelGeo", "oTe", MFnData::kMesh, &returnStatus);
	McheckErr(returnStatus, "ERROR creating lSystem output geometry attribute\n");
	numAttr.setStorable(false);


	returnStatus = addAttribute(MeshReduceNode::outputCurveCount);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");

	returnStatus = addAttribute(MeshReduceNode::inputFileName);
	McheckErr(returnStatus, "ERROR adding keys attribute\n");

	returnStatus = addAttribute(MeshReduceNode::outputCurves);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");

	returnStatus = addAttribute(MeshReduceNode::outputCtrlPointCount);
	McheckErr(returnStatus, "ERROR adding keys attribute\n");

	returnStatus = addAttribute(MeshReduceNode::outputJoints);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");

	returnStatus = addAttribute(MeshReduceNode::outputCurveIndices);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");


	returnStatus = addAttribute(MeshReduceNode::outputJointCount);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");

	returnStatus = addAttribute(MeshReduceNode::outputGeo);
	McheckErr(returnStatus, "ERROR adding output geo attribute\n");

	returnStatus = addAttribute(MeshReduceNode::outputTelGeo);
	McheckErr(returnStatus, "ERROR adding output geo attribute\n");

	attributeAffects(inputFileName, outputCurves);
	attributeAffects(inputFileName, outputJoints);
	attributeAffects(inputFileName, outputCurveIndices);
	attributeAffects(inputFileName, outputCurveCount);
	attributeAffects(inputFileName, outputJointCount);
	attributeAffects(inputFileName, outputCtrlPointCount);
	attributeAffects(inputFileName, outputGeo);
	attributeAffects(inputFileName, outputTelGeo);

	return MS::kSuccess;
}

MStatus MeshReduceNode::compute(const MPlug& plug, MDataBlock& data)
{
	MStatus status;

	if (plug == outputJoints || hasComputed)
	{
		MDataHandle geoData = data.inputValue(inputFileName, &status);
		McheckErr(status, "Error getting input mesh handle");

		MArrayDataHandle splineKeysData = data.outputArrayValue(outputCurves, &status);
		McheckErr(status, "Error getting output curves handle");

		MDataHandle curveIdxData = data.outputValue(outputCurveIndices, &status);
		McheckErr(status, "Error getting output curves indices handle");

		MArrayDataHandle jointsData = data.outputArrayValue(outputJoints, &status);
		McheckErr(status, "Error getting output joints handle");

		MDataHandle curveCountData = data.outputValue(outputCurveCount, &status);
		McheckErr(status, "Error getting output curves count handle");

		MDataHandle ctrlPointCountData = data.outputValue(outputCtrlPointCount, &status);
		McheckErr(status, "Error getting output curves count handle");

		MDataHandle jointCountData = data.outputValue(outputJointCount, &status);
		McheckErr(status, "Error getting output joints count handle");

		MDataHandle outGeoData = data.outputValue(outputGeo, &status);
		McheckErr(status, "Error getting output mesh value");

		MDataHandle outTelGeoData = data.outputValue(outputTelGeo, &status);
		McheckErr(status, "Error getting output mesh value");

		

		if (!hasComputed) 
		{
			halfEdgeMesh.isLoadObj = true;
			halfEdgeMesh.fileNameStd = geoData.asString().asChar();
			halfEdgeMesh.create();
			hasComputed = true;

			// Set outputMesh
			MPointArray pArray = MPointArray();
			MIntArray fConnect = MIntArray();
			MIntArray fCount = MIntArray();

			MFnMesh meshFS;

			float numVertices = halfEdgeMesh.vertices.size();
			for (int i = 0; i < numVertices; i++)
			{
				MPoint curPoint = MPoint(halfEdgeMesh.vertices.at(i)->pos.x,
					halfEdgeMesh.vertices.at(i)->pos.y,
					halfEdgeMesh.vertices.at(i)->pos.z,
					1.0);

				pArray.append(curPoint);
			}

			float numFaces = halfEdgeMesh.faces.size();
			for (int i = 0; i < numFaces; i++)
			{
				HalfEdge* he = halfEdgeMesh.faces.at(i)->edge;
				HalfEdge* first = he;

				int vertCount = 0;

				do
				{
					fConnect.append(he->nextVertex->id);
					vertCount++;
					he = he->nextEdge;
				} while (he != first);

				fCount.append(vertCount);
			}

			int a = fCount.length();
			int b = fConnect.length();

			MFnMeshData dataCreator;
			MObject newOutputData = dataCreator.create(&status);
			McheckErr(status, "ERROR creating outputData");

			MObject newMesh = meshFS.create(numVertices, numFaces,
				pArray, fCount, fConnect,
				newOutputData, &status);

			outGeoData.set(newOutputData);

			InitQuadrics();
			InitPQ();
			CollapseAll();

			int curveCount = outSplineKeys.size();
			int jointCount = outPoints.size();

			int pointCount = 0;
			for (int j = 0; j < outPoints.size(); j++) 
			{
				if ((!(outPoints.at(j)->pos.x - 0.0f) < FLT_EPSILON)
					&& (!(outPoints.at(j)->pos.y - 0.0f) < FLT_EPSILON)
					&& (!(outPoints.at(j)->pos.z - 0.0f) < FLT_EPSILON)) 
				{
					pointCount++;
				}
			}

			cubes = std::vector<std::vector<unique_ptr<Cube>>>(pointCount);

			for (int i = 0; i < curveCount; i++)
			{
				CubeArray ctrlGen = CubeArray();
				std::unique_ptr<Curve> curDiscreteCurve = std::make_unique<Curve>();
				for (int j = 0; j < outSplineKeys.at(i).size(); j++)
				{
					// Generate spline ctrl points
					glm::vec4 cKeyPos = glm::vec4(outSplineKeys.at(i).at(j), 1.0f);
					ctrlGen.keys.push_back(cKeyPos);
				}

				bool isHeadJoint = false;
				bool isTailJoint = false;
				int headJointIdx = -1;
				int tailJointIdx = -1;
				int pointCount = 0;
				for (int j = 0; j < outPoints.size(); j++) 
				{
					
					if ((!(outPoints.at(j)->pos.x - 0.0f) < FLT_EPSILON)
						&& (!(outPoints.at(j)->pos.y - 0.0f) < FLT_EPSILON)
						&& (!(outPoints.at(j)->pos.z - 0.0f) < FLT_EPSILON))
					{

						if (glm::length(outPoints.at(j)->pos - outSplineKeys.at(i).at(0)) < FLT_EPSILON) 
						{
							outPoints.at(j)->count++;
							isHeadJoint = true;
							headJointIdx = pointCount;
						}
						else if (glm::length(outPoints.at(j)->pos 
							     - outSplineKeys.at(i).at(outSplineKeys.at(i).size() - 1)) < FLT_EPSILON) 
						{
							outPoints.at(j)->count++;
							isTailJoint = true;
							tailJointIdx = pointCount;
						}
						pointCount++;
					}
				}

				if (isHeadJoint && !isTailJoint)
				{
					unique_ptr<Cube> headCube = make_unique<Cube>(glm::mat4(0.0f), Cube::JUNCTURE);
					unique_ptr<Cube> tailCube = make_unique<Cube>(glm::mat4(0.0f), Cube::GENERATOR);
					headCube->transform = glm::transpose(glm::translate(glm::mat4(), outSplineKeys.at(i).at(0)));
					tailCube->transform = glm::transpose(glm::translate(glm::mat4(),
						outSplineKeys.at(i).at(outSplineKeys.at(i).size() - 1)));
					curDiscreteCurve->parentCube = headCube.get();
					curDiscreteCurve->childCube = tailCube.get();
					headCube->parentCurve = curDiscreteCurve.get();
					tailCube->parentCurve = curDiscreteCurve.get();
					tailCube->rootCube = headCube.get();
					tailCube->parentJuncture = headCube.get();
					headCube->junctureId = headJointIdx;
					cubes.at(headJointIdx).push_back(std::move(headCube));
					otherCubes.push_back(std::move(tailCube));
				}
				else if (!isHeadJoint && isTailJoint)
				{
					unique_ptr<Cube> headCube = make_unique<Cube>(glm::mat4(0.0f), Cube::JUNCTURE);
					unique_ptr<Cube> tailCube = make_unique<Cube>(glm::mat4(0.0f), Cube::GENERATOR);
					tailCube->transform = glm::transpose(glm::translate(glm::mat4(), outSplineKeys.at(i).at(0)));
					headCube->transform = glm::transpose(glm::translate(glm::mat4(),
						outSplineKeys.at(i).at(outSplineKeys.at(i).size() - 1)));
					headCube->parentCurve = curDiscreteCurve.get();
					tailCube->parentCurve = curDiscreteCurve.get();
					curDiscreteCurve->parentCube = headCube.get();
					curDiscreteCurve->childCube = tailCube.get();
					headCube->rootCube = tailCube.get();
					headCube->parentJuncture = tailCube.get();
					tailCube->junctureId = tailJointIdx;
					cubes.at(tailJointIdx).push_back(std::move(tailCube));
					otherCubes.push_back(std::move(headCube));
				}
				else if (isHeadJoint && isTailJoint)
				{
					unique_ptr<Cube> headCube = make_unique<Cube>(glm::mat4(0.0f), Cube::JUNCTURE);
					unique_ptr<Cube> tailCube = make_unique<Cube>(glm::mat4(0.0f), Cube::JUNCTURE);
					tailCube->transform = glm::transpose(glm::translate(glm::mat4(), outSplineKeys.at(i).at(0)));
					headCube->transform = glm::transpose(glm::translate(glm::mat4(),
						outSplineKeys.at(i).at(outSplineKeys.at(i).size() - 1)));
					headCube->parentCurve = curDiscreteCurve.get();
					tailCube->parentCurve = curDiscreteCurve.get();
					curDiscreteCurve->parentCube = headCube.get();
					curDiscreteCurve->childCube = tailCube.get();
					headCube->junctureId = headJointIdx;
					tailCube->junctureId = tailJointIdx;
					headCube->neighbor = tailCube.get();
					tailCube->neighbor = headCube.get();
					cubes.at(headJointIdx).push_back(std::move(headCube));
					cubes.at(tailJointIdx).push_back(std::move(tailCube));
				}

				ctrlGen.computeCtrlPoints();
				ctrlGen.interpolate();
				curDiscreteCurve->points = new std::vector<glm::vec4>();
				for (int j = 0; j < ctrlGen.curve.size(); j++)
				{
					curDiscreteCurve->points->push_back(ctrlGen.curve.at(j));
					
				}
				curDiscreteCurve->reAssignPoints();
				curDiscreteCurve->discretilize();
				curDiscreteCurve->makeImpulseCurve();

				curDiscreteCurve->makeTelescope();
				delete(curDiscreteCurve->points);
				discreteCurves.push_back(std::move(curDiscreteCurve));				
			}

			for (int i = 0; i < cubes.size(); i++) 
			{
				for (int j = 0; j < cubes.at(i).size(); j++) 
				{
					if (cubes.at(i).at(j)->neighbor != nullptr) 
					{
						int neighborId = cubes.at(i).at(j)->neighbor->junctureId;
						int curId = cubes.at(i).at(j)->junctureId;
						int curNodeCount = outPoints.at(curId)->count;
						int neighborCount = outPoints.at(neighborId)->count;
						if (curNodeCount >= neighborCount)
						{
							if (!cubes.at(neighborId).at(0)->processed) 
							{
								for (int k = 0; k < cubes.at(neighborId).size(); k++)
								{
									cubes.at(neighborId).at(k)->parentJuncture = cubes.at(i).at(j).get();
									cubes.at(neighborId).at(k)->rootCube = cubes.at(i).at(j).get();
									cubes.at(neighborId).at(k)->processed = true;
								}
							}
							
							if (!cubes.at(i).at(j)->parentCurve->processed) 
							{
								cubes.at(i).at(j)->parentCurve->parentCube = cubes.at(i).at(j).get();
								cubes.at(i).at(j)->parentCurve->childCube = cubes.at(i).at(j)->neighbor;
								cubes.at(i).at(j)->parentCurve->processed = true;
							}
						}
						else if (curNodeCount < neighborCount) 
						{
							if (!cubes.at(i).at(0)->processed)
							{
								for (int k = 0; k < cubes.at(i).size(); k++)
								{
									cubes.at(i).at(k)->parentJuncture = cubes.at(i).at(j)->neighbor;
									cubes.at(i).at(k)->rootCube = cubes.at(i).at(j)->neighbor;
									cubes.at(i).at(k)->processed = true;
								}
							}

							if (!cubes.at(i).at(j)->parentCurve->processed)
							{
								cubes.at(i).at(j)->parentCurve->childCube = cubes.at(i).at(j).get();
								cubes.at(i).at(j)->parentCurve->parentCube = cubes.at(i).at(j)->neighbor;
								cubes.at(i).at(j)->parentCurve->processed = true;
							}
						}
					}
					
				}
				for (int i = 0; i < cubes.size(); i++)
				{

					if (cubes.at(i).at(0)->parentJuncture == nullptr) 
					{
						for (int j = 0; j < cubes.at(i).size(); j++) 
						{
							cubes.at(i).at(j)->parentJuncture = cubes.at(i).at(j).get();
						}
					}
				}
			}
		}
		

		// Output Data
		int curveCount = outSplineKeys.size();
		int jointCount = outPoints.size();

		curveCountData.setInt(curveCount);
		jointCountData.setInt(jointCount);

		int splineCount = 0;

		MIntArray curveIdx;


		// Set curveData
		// Compute control points for spline key
		
		std::vector<std::vector<glm::vec3>> ctrlPoints;
		for (int i = 0; i < curveCount; i++) 
		{
			CubeArray ctrlGen = CubeArray();
			for (int j = 0; j < outSplineKeys.at(i).size(); j++) 
			{
				// Generate spline ctrl points
				glm::vec4 cKeyPos = glm::vec4(outSplineKeys.at(i).at(j), 1.0f);
				ctrlGen.keys.push_back(cKeyPos);
			}
			ctrlGen.computeCtrlPoints();
			ctrlGen.interpolate();
			for (int j = 0; j < ctrlGen.curve.size(); j++) 
			{
				MDataHandle currentPoint = splineKeysData.outputValue();
				double4& point = currentPoint.asDouble4();
				point[0] = ctrlGen.curve.at(j).x;
				point[1] = ctrlGen.curve.at(j).y;
				point[2] = ctrlGen.curve.at(j).z;
				if (j == 0) 
				{
					point[3] = -1;
				}
				else 
				{
					point[3] = 1;
				}
				splineCount++;
				splineKeysData.next();
			}
		}
		
		ctrlPointCountData.setInt(splineCount);

		// Set jointData
		for (int i = 0; i < jointCount; i++)
		{
			MDataHandle currentJointData = jointsData.outputValue();
			double3& joint = currentJointData.asDouble3();
			joint[0] = outPoints.at(i)->pos.x;
			joint[1] = outPoints.at(i)->pos.y;
			joint[2] = outPoints.at(i)->pos.z;

			jointsData.next();
		}

		if (enterOnce) 
		{
			enterOnce = false;
			MFnMesh meshFS;

			float numVertices = 0.0f;
			float numFaces = 0.0f;
			MPointArray pArray;
			MIntArray fConnect;
			MIntArray fCount;
			int faceCounts;

			for (unique_ptr<Curve>& dc : discreteCurves)
			{
				dc->updateSegmentTransforms();
				for (unique_ptr<Shell>& pShell : dc->shells)
				{
					numFaces += pShell->ib.size() / 3.0f;
					for (MPoint p : pShell->points)
					{
						glm::vec4 curPoint = glm::vec4(p.x, p.y, p.z, 1.0f);
						curPoint = pShell->animatedTransform * curPoint;
						pArray.append(MPoint(curPoint.x, curPoint.y, curPoint.z));
					}
					for (int connection : pShell->faceConnects)
					{
						fConnect.append(connection + numVertices);
					}
					for (int count : pShell->faceCounts)
					{
						fCount.append(count);
					}

					faceCounts += pShell->faceCounts.length();
					numVertices += pShell->points.length();
				}
			}

			MFnMeshData dataCreator;
			MObject newOutputData = dataCreator.create(&status);
			McheckErr(status, "ERROR creating outputData");

			MObject newMesh = meshFS.create(numVertices, numFaces, pArray,
				fCount, fConnect,
				newOutputData, &status);

			outTelGeoData.set(newOutputData);
		}
		
	
		
		
		ctrlPointCountData.setClean();
		outGeoData.setClean();
		outTelGeoData.setClean();
		geoData.setClean();
		curveCountData.setClean();
		jointCountData.setClean();
		splineKeysData.setClean();
		jointsData.setClean();
		curveIdxData.setClean();
		data.setClean(plug);
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}

void MeshReduceNode::InitQuadrics() 
{
	for (int i = 0; i < halfEdgeMesh.faces.size(); i++)
	{
		Face* curFace = halfEdgeMesh.faces.at(i).get();

		if (isBoundary(curFace))
		{
			continue;
		}

		Vertex* faceVertex = curFace->edge->nextVertex;



		// Get a point on current face and normal vector
		glm::vec3 faceVertexPos = glm::vec3(faceVertex->pos);

		glm::vec3 faceVert1 = glm::vec3(curFace->edge->nextEdge->nextVertex->pos);
		glm::vec3 faceVert2 = glm::vec3(curFace->edge->nextEdge->nextEdge->nextVertex->pos);

		glm::vec3 vec1 = faceVert2 - faceVert1;
		glm::vec3 vec2 = faceVertexPos - faceVert1;

		glm::vec3 faceNormal = -glm::vec3(glm::normalize(glm::cross(vec1, vec2)));
		curFace->normal = glm::vec4(faceNormal, 0.0f);
		float d = -glm::dot(faceNormal, faceVertexPos);

		glm::vec4 v = glm::vec4(faceNormal.x, faceNormal.y, faceNormal.z, d);

		faceQuadric.push_back(glm::outerProduct(v, v));
	}

	uf = UnionFind();
	uf.isVertex = true;
	for (int i = 0; i < halfEdgeMesh.vertices.size(); i++) 
	{
		std::unique_ptr<VertexNode> curUnionVert = std::make_unique<VertexNode>();
		Vertex* curVert = halfEdgeMesh.vertices.at(i).get();

		curUnionVert->pos = glm::vec3(curVert->pos);
		curUnionVert->contaningUf = &uf;


		HalfEdge* he = curVert->edge;
		HalfEdge* first = he;
		
		glm::mat4 result = glm::mat4(0.0f);

		do 
		{
			Face* f = he->face;
			result += faceQuadric.at(f->id);
			he = he->symEdge->nextEdge;
		} while (he != first);

		vertQuadric.push_back(result);

		curUnionVert->quadric = result;
		curUnionVert->mergedPoints = std::unordered_set<int>();
		curUnionVert->mergedPoints.insert(i);
		curUnionVert->index = i;

		uf.nodes.push_back(std::move(curUnionVert));
	}

	for (int i = 0; i < halfEdgeMesh.halfEdges.size(); i++)
	{
		glm::mat4 Qsum = glm::mat4(0.0f);
		int v1 = halfEdgeMesh.halfEdges.at(i)->nextVertex->id;
		int v2 = halfEdgeMesh.halfEdges.at(i)->symEdge->nextVertex->id;
		Qsum += vertQuadric.at(v1);
		Qsum += vertQuadric.at(v2);

		glm::vec3 opt = optimalPoint(Qsum, glm::vec3(halfEdgeMesh.vertices.at(v1)->pos), glm::vec3(halfEdgeMesh.vertices.at(v2)->pos));
		edgeCosts.push_back(errorAtPoint(Qsum, opt));
	}
}

void MeshReduceNode::InitPQ()
{
	PQ = priority_queue<FaceNode*, std::vector<FaceNode*>, cmp>();

	for (int i = 0; i < halfEdgeMesh.faces.size(); i++) 
	{
		Face* f = halfEdgeMesh.faces.at(i).get();
		HalfEdge* minEdge = optimalEdgeShortest(f);
		float cost = lengthCost(minEdge);
		std::unique_ptr<FaceNode> curFaceNode = std::make_unique<FaceNode>(f, cost);
		PQ.push(curFaceNode.get());
		faceNodes.push_back(std::move(curFaceNode));
	}
}

void MeshReduceNode::CollapseAll() 
{
	while (PQ.size() != 0) 
	{
		FaceNode* nextFace = PQ.top();
		PQ.pop();

		if (nextFace->face->collapsed) continue;

		Face* f = nextFace->face;

		HalfEdge* minEdge = optimalEdgeShortest(f);

		// Collapse Edge
		int v1 = minEdge->nextVertex->id; 
		int v2 = minEdge->symEdge->nextVertex->id;

		minEdge->face->collapsed = true;
		minEdge->symEdge->face->collapsed = true;

		uf.Union(v1, v2);

		std::unordered_set<int> vPoints = dynamic_cast<VertexNode*>(uf.find(v1))->mergedPoints;

		for (auto iter = vPoints.begin(); iter != vPoints.end(); ++iter)
		{
			halfEdgeMesh.vertices.at((int)(*iter))->pos = glm::vec4(dynamic_cast<VertexNode*>(uf.find(v1))->pos, 1.0f);
		}
		
		// Update Priorities
		for (auto iter = vPoints.begin(); iter != vPoints.end(); ++iter)
		{
			FixPQ(halfEdgeMesh.vertices.at((int)(*iter)).get());
		}
	}

	// Make skeleton
	// Adjecency lists
	adjacency = std::vector<std::unordered_set<int>*>(halfEdgeMesh.vertices.size());
	std::vector<bool> isCreated = std::vector<bool>(halfEdgeMesh.vertices.size());
	for (int i = 0; i < halfEdgeMesh.vertices.size(); i++) 
	{
		isCreated.at(i) = false;
	}

	for (int i = 0; i < halfEdgeMesh.vertices.size(); i++) 
	{
		VertexNode* root = dynamic_cast<VertexNode*>(uf.find(i));
		unordered_set<int> points = root->mergedPoints;

		if (!isCreated.at(root->index))
		{
			isCreated.at(root->index) = true;
			std::unique_ptr<std::unordered_set<int>> adj = std::make_unique<std::unordered_set<int>>();
			for (auto iter = points.begin(); iter != points.end(); ++iter) 
			{
				adjacency.at((int)(*iter)) = adj.get();
			}
			adjStore.push_back(std::move(adj));
		}

		for (auto iter = points.begin(); iter != points.end(); ++iter) 
		{
			Vertex* curVert = halfEdgeMesh.vertices.at((int)(*iter)).get();
			HalfEdge* startEdge = curVert->edge;
			HalfEdge* he = startEdge;

			do 
			{
				int neighborInd = he->nextVertex->id;
				int neighborRoot = dynamic_cast<VertexNode*>(uf.find(neighborInd))->index;

				if (points.find(neighborRoot) == points.end()) 
				{
					adjacency.at((int)(*iter))->insert(neighborRoot);
				}

				he = he->symEdge->nextEdge;

			} while (he != startEdge);
		}
	}
		
	// Find highest degree
	int maxDeg = 0;
	int maxVert = 0;
	for (int i = 0; i < halfEdgeMesh.vertices.size(); i++) 
	{
		int degree = adjacency.at(i)->size();
		if (degree > maxDeg) 
		{
			maxVert = i;
			maxDeg = degree;
		}
	}

	std::vector<bool> visited = std::vector<bool>(halfEdgeMesh.vertices.size());
	for(int i = 0; i < visited.size(); i++)
	{
		visited.at(i) = false;
	}

	std::stack<IntPair> dfsStack;

	IntPair firstEntry(maxVert, maxVert);

	dfsStack.push(firstEntry);

	std::list<int> curPoints;

	std::map<int, BulbPoints*> bulbDict;
	BulbPoints* firstBulb;

	while (dfsStack.size() > 0) 
	{
		int2 nextPair = { dfsStack.top().first, dfsStack.top().second };
		dfsStack.pop();
		int next = nextPair[0];

		if (visited[next]) continue;

		std::unordered_set<int> equivClass = dynamic_cast<VertexNode*>(uf.find(next))->mergedPoints;

		for(auto iter = equivClass.begin(); iter != equivClass.end(); ++iter)
		{
			visited.at((int)(*iter)) = true;
		}

		std::unordered_set<int>* neighbors = adjacency.at(next);
		for (auto iter = neighbors->begin(); iter != neighbors->end(); ++iter)
		{
			if(!visited.at(int(*iter)))
				dfsStack.push(IntPair(int(*iter), next));
		}

		if (curPoints.empty()) 
		{
			// Transform not added
			glm::vec3 pos = glm::vec3(halfEdgeMesh.vertices.at(next)->pos);
			std::unique_ptr<BulbPoints> firstBulbPtr = std::make_unique<BulbPoints>();
			firstBulbPtr->pos = pos;
			firstBulbPtr->radius = 0.5f;
			firstBulb = firstBulbPtr.get();
			outPoints.push_back(std::move(firstBulbPtr));
			

			bulbDict[next] = firstBulb;

			curPoints.push_back(next);
		}
		else if (neighbors->size() == 2) 
		{
			// Degree equals two
			curPoints.push_back(next);
		}
		else if (neighbors->size() == 1)
		{
			curPoints.push_back(next);

			std::vector<glm::vec3> spline;

			for (auto iter = curPoints.begin(); iter != curPoints.end(); iter++)
			{
				spline.push_back(glm::vec3(halfEdgeMesh.vertices.at(int(*iter))->pos));
			}

			outSplineKeys.push_back(spline);

			curPoints.clear();

			firstBulb->index = outSplineKeys.size() - 1;
			firstBulb->isFirst = true;

			if (dfsStack.size() > 0)
			{
				curPoints.push_back(dfsStack.top().second);
				std::unique_ptr<BulbPoints> bulb;
				if (bulbDict.find(dfsStack.top().second) != bulbDict.end())
				{
					bulb = std::make_unique<BulbPoints>();
					firstBulb = bulb.get();
					outPoints.push_back(std::move(bulb));
				}
				else
					firstBulb = nullptr;
			}
		}
		else 
		{
			// Transform not add
			glm::vec3 pos = glm::vec3(halfEdgeMesh.vertices.at(next)->pos);
			std::unique_ptr<BulbPoints> bulb = std::make_unique<BulbPoints>();
			bulb->pos = pos;
			bulb->radius = 0.5f;
			bulbDict[next] = bulb.get();

			curPoints.push_back(next);
			std::vector<glm::vec3> spline;

			for (auto iter = curPoints.begin(); iter != curPoints.end(); iter++)
			{
				spline.push_back(glm::vec3(halfEdgeMesh.vertices.at(int(*iter))->pos));
			}

			outSplineKeys.push_back(spline);
			bulb->index = outSplineKeys.size() - 1;
			firstBulb->index = outSplineKeys.size() - 1;
			firstBulb = bulb.get();
			outPoints.push_back(std::move(bulb));
			curPoints.clear();
			curPoints.push_back(next);
		}
	}
}

HalfEdge* MeshReduceNode::optimalEdge(Face* f) 
{
	HalfEdge* e1 = f->edge;
	HalfEdge* e2 = f->edge->nextEdge;
	HalfEdge* e3 = f->edge->nextEdge->nextEdge;

	HalfEdge* minEdge = e1;
	float minCost = edgeCosts.at(e1->id);
	float cost2 = edgeCosts.at(e2->id);
	float cost3 = edgeCosts.at(e3->id);

	if (cost2 < minCost)
	{
		minCost = cost2;
		minEdge = e2;
	}

	if (cost3 < minCost)
	{
		minCost = cost3;
		minEdge = e3;
	}
	return minEdge;
}

float MeshReduceNode::Cost(HalfEdge* he) 
{
	int v1 = he->nextVertex->id;
	int v2 = he->symEdge->nextVertex->id;

	glm::mat4 Qsum = glm::mat4(0.0f);
	Qsum += vertQuadric.at(v1);
	Qsum += vertQuadric.at(v2);

	glm::vec3 opt = optimalPoint(Qsum,
		glm::vec3(he->nextVertex->pos),
		glm::vec3(he->symEdge->nextVertex->pos));

	float cost = errorAtPoint(Qsum, opt);
	return cost;
}

bool MeshReduceNode::isBoundary(Face* f) 
{
	HalfEdge* he = f->edge;
	HalfEdge* first = he;
	do 
	{
		if (he->symEdge == nullptr)
			return true;
		he = he->nextEdge;
	} while (he != first);

	return false;
}

void MeshReduceNode::FixPQ(Vertex* v) 
{
	HalfEdge* start = v->edge;
	HalfEdge* he = start;

	// Update edge costs for affected (adjacent) edges
	do
	{
		edgeCosts[he->id] = edgeQuadric(he);
		edgeCosts[he->symEdge->id] = edgeCosts[he->id];

		he = he->symEdge->nextEdge;
	} while (he != start);

	start = v->edge;
	he = start;

	// Iterate over all adjacent edges
	do
	{
		Face* f = he->face;
		if (!f->collapsed)
		{
			FaceNode* updateNode = faceNodes.at(f->id).get();
			HalfEdge* minEdge = optimalEdgeShortest(f);
			updateNode->cost = lengthCost(minEdge);

			FaceNode* dummy = PQ.top();
			PQ.pop();
			PQ.push(dummy);
		}

		he = he->symEdge->nextEdge;
	} while (he != start);
}

float MeshReduceNode::edgeQuadric(HalfEdge* he) 
{
	int v1 = he->symEdge->nextVertex->id;
	int v2 = he->nextVertex->id;

	glm::mat4 Qsum = glm::mat4(0.0f);
	Qsum += dynamic_cast<VertexNode*>(uf.find(v1))->quadric;
	Qsum += dynamic_cast<VertexNode*>(uf.find(v2))->quadric;

	glm::vec3 opt = optimalPoint(Qsum, glm::vec3(he->nextVertex->pos), glm::vec3(he->symEdge->nextVertex->pos));
	return errorAtPoint(Qsum, opt);
}

HalfEdge* MeshReduceNode::optimalEdgeShortest(Face* f)
{
	HalfEdge* e1 = f->edge;
	HalfEdge* e2 = e1->nextEdge;
	HalfEdge* e3 = e2->nextEdge;

	HalfEdge* shortestHE = e1;
	float shortest = glm::length(e1->nextVertex->pos - e1->symEdge->nextVertex->pos);
	float e2len = glm::length(e2->nextVertex->pos - e2->symEdge->nextVertex->pos);
	float e3len = glm::length(e3->nextVertex->pos - e3->symEdge->nextVertex->pos);

	if (e2len < shortest)
	{
		shortest = e2len;
		shortestHE = e2;
	}
	if (e3len < shortest)
	{
		shortest = e3len;
		shortestHE = e3;
	}

	return shortestHE;
}

float MeshReduceNode::lengthCost(HalfEdge* he) 
{
	return glm::length(he->nextVertex->pos - he->symEdge->nextVertex->pos);
}