#include "CurveGenNode.h"
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

MTypeId CurveGenNode::id(0x80001);
MObject CurveGenNode::inKeys;
MObject CurveGenNode::outCtrlPoint;

void* CurveGenNode::creator() 
{
	return new CurveGenNode;
}

MStatus CurveGenNode::initialize() 
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
	
	CurveGenNode::inKeys = numAttr.create("splineKeys", "key", MFnNumericData::k3Double, 0.0f, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve keys attribute");

	numAttr.setArray(true);

	numAttr.setKeyable(true);
	numAttr.setStorable(true);
	numAttr.setReadable(true);
	numAttr.setWritable(true);

	CurveGenNode::outCtrlPoint = numAttr.create("ctrlPoints", "cp", MFnNumericData::k3Double, 0.0f, &returnStatus);
	McheckErr(returnStatus, "ERROR creating curve control points attribute");

	numAttr.setArray(true);
	numAttr.setStorable(false);

	returnStatus = addAttribute(CurveGenNode::inKeys);
	McheckErr(returnStatus, "ERROR adding keys attribute\n");

	returnStatus = addAttribute(CurveGenNode::outCtrlPoint);
	McheckErr(returnStatus, "ERROR adding control points attribute\n");

	attributeAffects(inKeys, outCtrlPoint);

	return MS::kSuccess;
}

MStatus CurveGenNode::compute(const MPlug& plug, MDataBlock& data) 
{
	MStatus status;

	if (plug == outCtrlPoint)
	{
		MArrayDataHandle keysData = data.inputArrayValue(inKeys, &status);
		McheckErr(status, "Error getting time handle");
		int keysCount = keysData.elementCount();

		MArrayDataHandle ctrlData = data.outputArrayValue(outCtrlPoint, &status);
		McheckErr(status, "Error getting output geo handle");

		// Set curve keys to compute ctrl point
		CubeArray ctrlGen = CubeArray();

		for (int i = 0; i < keysCount; i++) 
		{
			MDataHandle currentKey = keysData.inputValue();
			double3& keyPos = currentKey.asDouble3();
			glm::vec4 cKeyPos = glm::vec4(keyPos[0], keyPos[1], keyPos[2], 1.0f);
			ctrlGen.keys.push_back(cKeyPos);

			if (!keysData.next()) break;
		}

		ctrlGen.computeCtrlPoints();
		ctrlGen.interpolate();

		for (int i = 0; i < ctrlGen.curve.size(); i++)
		{
			MDataHandle currentPoint = ctrlData.outputValue();
			double3& point = currentPoint.asDouble3();
			point[0] = ctrlGen.curve[i].x;
			point[1] = ctrlGen.curve[i].y;
			point[2] = ctrlGen.curve[i].z;

			ctrlData.next();
		}

		keysData.setClean();
		ctrlData.setClean();
		data.setClean(plug);
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}