#pragma once
#include <maya/MPxNode.h>
#include <maya/MStatus.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MDataBlock.h>
#include <maya/MObject.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>

#define McheckErr(stat,msg)                     \
        if ( MS::kSuccess != stat ) {   \
                cerr << msg;                            \
                return MS::kFailure;            \
        }

class CurveGenNode : public MPxNode
{
public:
	CurveGenNode() {};
	virtual ~CurveGenNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static void* creator();
	static MStatus initialize();

	static MTypeId id;

	// Telescope Attributes
	//static MObject impulseCount;
	static MObject inKeys;
	static MObject outCtrlPoint;
};