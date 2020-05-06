#include <maya/MPxCommand.h>
#include <maya/MFnPlugin.h>
#include <maya/MIOStream.h>
#include <maya/MString.h>
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSimple.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MDGModifier.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MStringArray.h>
#include "CurveGenNode.h"
#include "curveOptNode.h"
#include "meshReduceNode.h"
#include <list>

#include "LSystemCmd.h"

MStatus initializePlugin( MObject obj )
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj, "MyPlugin", "1.0", "Any");


    // Register Node
    status = plugin.registerNode("CurveGenNode", CurveGenNode::id, CurveGenNode::creator, CurveGenNode::initialize);
    status = plugin.registerNode("CurveOptNode", CurveOptNode::id, CurveOptNode::creator, CurveOptNode::initialize);
    status = plugin.registerNode("MeshReduceNode", MeshReduceNode::id, MeshReduceNode::creator, MeshReduceNode::initialize);
    if (!status)
    {
        status.perror("registerNode");
        return status;
    }

    return status;
}

MStatus uninitializePlugin( MObject obj)
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj );


    status = plugin.deregisterNode(CurveGenNode::id);
    status = plugin.deregisterNode(CurveOptNode::id);
    status = plugin.deregisterNode(MeshReduceNode::id);
    if (!status)
    {
        status.perror("deregisterNode");
        return status;
    }

    return status;
}


