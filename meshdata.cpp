#include "meshdata.h"
#include "global.h"

int Vertex::lastId = 0;

Vertex::Vertex():
    isSharped(false),
    sharpCount(0),
    sharpness(1.0f)
{
    id = lastId;
    lastId++;
}

Vertex::Vertex(int _id):
    id(_id) {}

Vertex::Vertex(glm::vec4& _pos, std::unique_ptr<HalfEdge> _edge, int _id):
    pos(_pos), edge(_edge.get()), id(_id) {}


int HalfEdge::lastId = 0;

HalfEdge::HalfEdge():
    isSharped(false),
    sharpness(1.0f)
{
    id = lastId;
    lastId++;
}

HalfEdge::HalfEdge(int _id):
    id(_id) {}

HalfEdge::HalfEdge(std::unique_ptr<HalfEdge> _nextEdge,
                   std::unique_ptr<HalfEdge> _symEdge,
                   std::unique_ptr<Face> _face,
                   std::unique_ptr<Vertex> _nextVertex, int _id):
    nextEdge(_nextEdge.get()), symEdge(_symEdge.get()),
    face(_face.get()), nextVertex(_nextVertex.get()), id(_id) {}

int Face::lastId = 0;

Face::Face():
    isMid(false),
    isSharped(false),
    sharpness(1.0f),
    collapsed(false)
{
    id = lastId;
    lastId++;
}

Face::Face(int _id):
    id(_id) {}

Face::Face(glm::vec4& _color, std::unique_ptr<HalfEdge> _edge, int _id):
    color(_color), edge(_edge.get()), id(_id) {}

UnionFindNode* UnionFind::find(UnionFindNode* node) 
{
    if (node->parent != node) 
    {
        node->parent = find(node->parent);
    }
    return node->parent;
}

UnionFindNode* UnionFind::find(int x) 
{
    return find(nodes.at(x).get());
}

void UnionFind::Union(UnionFindNode* x, UnionFindNode* y) 
{
    UnionFindNode* xRoot = find(x);
    UnionFindNode* yRoot = find(y);

    if (xRoot == yRoot) return;

    // Union-by-rank: make the smaller set the child of the larger set.
    if (xRoot->rank < yRoot->rank)
    {
        yRoot->Merge(xRoot);
        xRoot->parent = yRoot;
    }
    else if (xRoot->rank > yRoot->rank)
    {
        xRoot->Merge(yRoot);
        yRoot->parent = xRoot;
    }
    else
    {
        xRoot->Merge(yRoot);
        // If same size, break ties either way and increase rank.
        yRoot->parent = xRoot;
        xRoot->rank++;
    }
}

void UnionFind::Union(int x, int y)
{
    Union(nodes.at(x).get(), nodes.at(y).get());
}

void VertexNode::Merge(UnionFindNode* child) 
{
    VertexNode* childNode = dynamic_cast<VertexNode*>(child);

    if (child->parent != child) 
        std::cout << "No parent for this node!" << std::endl;

    quadric += childNode->quadric;

    float n1 = mergedPoints.size();
    float n2 = childNode->mergedPoints.size();

    pos = (n1 * pos + n2 * childNode->pos) / (n1 + n2);

    currentCost = errorAtPoint(quadric, pos);

    for (auto iter = childNode->mergedPoints.begin();
        iter != childNode->mergedPoints.end();
        ++iter) 
    {
        mergedPoints.insert((int)(*iter));
    }
}