#include "curveOptNode.h"
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
#include "cubeArray.h"
#include "curve.h"
#include "cube.h"

MTypeId CurveOptNode::id(0x80002);
MObject CurveOptNode::inPoints;
MObject CurveOptNode::outOptPoints;
MObject CurveOptNode::outPointsCount;
MObject CurveOptNode::isAdding;
MObject CurveOptNode::outputGeo;
MObject CurveOptNode::time;
MObject CurveOptNode::isExtending;
MObject CurveOptNode::isAdjusting;


void* CurveOptNode::creator()
{
	return new CurveOptNode;
}

MStatus CurveOptNode::initialize()
{
	MFnNumericAttribute numAttr;
	MFnUnitAttribute unitAttr;
	MFnTypedAttribute typeAttr;
	MStatus returnStatus;

	cout << "InitializeNode!" << endl;

	MVectorArray defaultPoints;
	MFnVectorArrayData defaultArray;
	MObject defaultAttr;

	defaultPoints.clear();
	defaultAttr = defaultArray.create(defaultPoints);

	CurveOptNode::inPoints = numAttr.create("curvePoints", "cps", MFnNumericData::k3Double, 0.0f, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve keys attribute");

	numAttr.setArray(true);

	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	CurveOptNode::time = unitAttr.create("time", "tm", MFnUnitAttribute::kTime, 0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating lSystem time attribute\n");
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	CurveOptNode::isAdding = numAttr.create("isAdding", "ad", MFnNumericData::kBoolean, false, &returnStatus);
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	CurveOptNode::isExtending = numAttr.create("isExtending", "ise", MFnNumericData::kBoolean, true, &returnStatus);
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);


	CurveOptNode::isAdjusting = numAttr.create("isAdjusting", "iaj", MFnNumericData::kBoolean, true, &returnStatus);
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	CurveOptNode::outPointsCount = numAttr.create("pointsCount", "pc", MFnNumericData::kInt, 0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating points count attribute");
	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	CurveOptNode::outputGeo = typeAttr.create("outGeo", "oGe", MFnData::kMesh, &returnStatus);
	McheckErr(returnStatus, "ERROR creating lSystem output geometry attribute\n");
	numAttr.setStorable(false);

	CurveOptNode::outOptPoints = numAttr.create("optPoints", "op", MFnNumericData::k3Double, 0.0f, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setArray(true);
	numAttr.setStorable(false);

	returnStatus = addAttribute(CurveOptNode::inPoints);
	McheckErr(returnStatus, "ERROR adding keys attribute\n");

	returnStatus = addAttribute(CurveOptNode::time);
	McheckErr(returnStatus, "ERROR adding keys attribute\n");

	returnStatus = addAttribute(CurveOptNode::outPointsCount);
	McheckErr(returnStatus, "ERROR adding points count attribute\n");

	returnStatus = addAttribute(CurveOptNode::outOptPoints);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");

	returnStatus = addAttribute(CurveOptNode::isAdding);
	McheckErr(returnStatus, "ERROR adding isAdding attribute\n");

	returnStatus = addAttribute(CurveOptNode::isExtending);
	McheckErr(returnStatus, "ERROR adding isAdding attribute\n");

	returnStatus = addAttribute(CurveOptNode::isAdjusting);
	McheckErr(returnStatus, "ERROR adding isAdding attribute\n");

	returnStatus = addAttribute(CurveOptNode::outputGeo);
	McheckErr(returnStatus, "ERROR adding output geo attribute\n");

	attributeAffects(inPoints, outPointsCount);
	attributeAffects(inPoints, outOptPoints);
	attributeAffects(inPoints, outputGeo);
	attributeAffects(time, outputGeo);
	attributeAffects(isExtending, outputGeo);
	attributeAffects(isAdjusting, outputGeo);


	return MS::kSuccess;
}

MStatus CurveOptNode::compute(const MPlug& plug, MDataBlock& data)
{
	MStatus status;

	MDataHandle addingBool = data.inputValue(isAdding, &status);
	McheckErr(status, "ERROR getting isAdding handle");
	bool addBool = addingBool.asBool();
	
	if (plug == outputGeo)
	{
		MArrayDataHandle keysData = data.inputArrayValue(inPoints, &status);
		McheckErr(status, "Error getting time handle");
		int keysCount = keysData.elementCount();

		MDataHandle countData = data.outputValue(outPointsCount, &status);
		McheckErr(status, "Error getting output count handle");

		MArrayDataHandle ctrlData = data.outputArrayValue(outOptPoints, &status);
		McheckErr(status, "Error getting output geo handle");

		MDataHandle geoData = data.outputValue(outputGeo, &status);
		McheckErr(status, "Error getting output geo handle");

		MDataHandle timeData = data.inputValue(time, &status);
		McheckErr(status, "Error getting timer handle");

		MDataHandle extendingData = data.inputValue(isExtending, &status);
		McheckErr(status, "Error getting extendingData");

		MDataHandle adjustingData = data.inputValue(isExtending, &status);
		McheckErr(status, "Error getting extendingData");


		bool isExtend = extendingData.asBool();
		MFnMesh meshFS;

		bool isAdjusting = adjustingData.asBool();
		
		if (isAdjusting) 
		{
			discreteCurve.computed = false;
		}

		if (!discreteCurve.computed) 
		{
			// Set curve keys to compute ctrl point
			discreteCurve = Curve(nullptr, nullptr);
			discreteCurve.computed = true;
			discreteCurve.points = new std::vector<glm::vec4>();

			
			firstCube = Cube(glm::mat4(0.0f), 5);
			lastCube = Cube(glm::mat4(0.0f), 4);
			firstCube.parentCurve = nullptr;
			lastCube.parentCurve = nullptr;
			firstCube.parentJuncture = &firstCube;
			lastCube.rootCube = &firstCube;
			firstCube.rootCube = &firstCube;
			lastCube.parentJuncture = nullptr;


			for (int i = 0; i < keysCount; i++)
			{
				MDataHandle currentKey = keysData.inputValue();
				double3& keyPos = currentKey.asDouble3();
				glm::vec4 cKeyPos = glm::vec4(keyPos[0], keyPos[1], keyPos[2], 1.0f);
				if (i == 0)
				{
					glm::mat4 cubeMat = glm::translate(glm::mat4(0.0f), glm::vec3(cKeyPos));
					firstCube.transform = glm::transpose(cubeMat);
					discreteCurve.parentCube = &firstCube;
				}
				else if (i == keysCount - 1)
				{
					glm::mat4 cubeMat = glm::translate(glm::mat4(0.0f), glm::vec3(cKeyPos));
					lastCube.transform = glm::transpose(cubeMat);
					discreteCurve.childCube = &lastCube;
				}
				discreteCurve.points->push_back(cKeyPos);

				if (!keysData.next()) break;
			}

			discreteCurve.reAssignPoints();
			discreteCurve.discretilize();
			discreteCurve.makeImpulseCurve();

			int torsionPointSize = discreteCurve.torsionImpulsePoints.size();
			countData.set(torsionPointSize);

			discreteCurve.makeTelescope();
		}
		
		
		if (!isAdjusting) 
		{
			if (isExtend && (discreteCurve.extensionState == Curve::RETRACTED))
			{
				discreteCurve.extensionState = Curve::EXTENDING;
			}
			else if (discreteCurve.extensionState == Curve::EXTENDED)
			{
				if (!isExtend)
				{
					discreteCurve.extensionState = Curve::RETRACTING;
				}

			}

			if (discreteCurve.extensionState == Curve::EXTENDING)
			{
				discreteCurve.extensionExtent += 6.f / Curve::numImpulses / 1000.0f;
				if (discreteCurve.extensionExtent > 1.f)
				{
					discreteCurve.extensionState = Curve::EXTENDED;
					discreteCurve.extensionExtent = 1.f;
				}
			}
			else if (discreteCurve.extensionState == Curve::RETRACTING)
			{
				discreteCurve.extensionExtent -= 6.f / Curve::numImpulses / 1000.0f;
				if (discreteCurve.extensionExtent < 0.f)
				{
					discreteCurve.extensionState = Curve::RETRACTED;
					discreteCurve.extensionExtent = 0.f;
				}
			}
		}
		
		currentExtent = discreteCurve.extensionExtent;
		currentState = discreteCurve.extensionState;

		discreteCurve.updateSegmentTransforms();
		float numVertices = 0.0f;
		float numFaces = 0.0f;
		MPointArray pArray;
		MIntArray fConnect;
		MIntArray fCount;
		int faceCounts;
		for (unique_ptr<Shell>& pShell : discreteCurve.shells) 
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


		

		MFnMeshData dataCreator;
		MObject newOutputData = dataCreator.create(&status);
		McheckErr(status, "ERROR creating outputData");

		MObject newMesh = meshFS.create(numVertices, numFaces, pArray,
			                            fCount, fConnect, 
			                            newOutputData, &status);

		geoData.set(newOutputData);
		
		countData.setClean();
		keysData.setClean();
		ctrlData.setClean();
		data.setClean(plug);
		//delete(discreteCurve.points);
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}