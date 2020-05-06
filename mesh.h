#ifndef MESH
#define MESH
#include <vector>
#include "meshdata.h"
#include <glm/glm.hpp>

class Mesh
{
public:
    std::vector<std::unique_ptr<Face>> faces;
    std::vector<std::unique_ptr<HalfEdge>> halfEdges;
    std::vector<std::unique_ptr<Vertex>> vertices;
    HalfEdge* heRoot;
    glm::mat4 curModel;
    bool isLoadObj;
    std::string fileNameStd;
    int faceVerticesCount;
    int selectedId;

    int count;

    bool weightBound;
    bool weightIdxBound;

public:
    // Constructor
    Mesh();

    // Construct a mesh from file
    void constructFileMesh(std::string& fileName);

    // Destructor
    ~Mesh();
    // Create function
    void create();

    // Update data
    void update();
};

#endif // MESH

