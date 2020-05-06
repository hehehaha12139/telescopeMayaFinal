#pragma once
#include <maya/MPxNode.h>
#include <maya/MStatus.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MDataBlock.h>
#include <maya/MObject.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include "curve.h"

#define McheckErr(stat,msg)                     \
        if ( MS::kSuccess != stat ) {   \
                cerr << msg;                            \
                return MS::kFailure;            \
        }

class CurveOptNode : public MPxNode
{
public:
	CurveOptNode() {};
	virtual ~CurveOptNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static void* creator();
	static MStatus initialize();
	Curve discreteCurve;
	Cube firstCube;
	Cube lastCube;
	float currentExtent;
	int currentState;
	

	static MTypeId id;

	// Telescope Attributes
	//static MObject impulseCount;
	static MObject inPoints;
	static MObject outOptPoints;
	static MObject outPointsCount;
	static MObject isAdding;
	static MObject outputGeo;
	static MObject time;
	static MObject isExtending;
	static MObject isAdjusting;
};