#include "mesh.h"
#include "global.h"
#include <cmath>
#include <fStream>
#include <map>
#include <unordered_map>
#include <sstream>
#include <random>
#include <ctime>
#include <array>

Mesh::Mesh():
    isLoadObj(false), selectedId(-1) {}

void Mesh::create()
{
    // Construct half-edge data structure
    if(!isLoadObj)
    {
        Vertex::lastId = 0;
        HalfEdge::lastId = 0;
        Face::lastId = 0;
    }
    else
    {
        Vertex::lastId = 0;
        HalfEdge::lastId = 0;
        Face::lastId = 0;
        constructFileMesh(fileNameStd);
    }

    // Set data vector
    std::vector<glm::vec4> pos;
    std::vector<glm::vec4> col;
    std::vector<glm::vec4> nor;
    std::vector<unsigned int> idx;

    count = 0;

    for(int i = 0; i <  faces.size(); i++)
    {
        HalfEdge* first = faces[i]->edge;
        unsigned int prePosCount = pos.size();
        do {
            Vertex* curVertex = first->nextVertex;

            Vertex* lastVert = first->symEdge->nextVertex;
            Vertex* nextVert = first->nextEdge->nextVertex;
            glm::vec3 edge1 = glm::vec3(curVertex->pos - lastVert->pos);
            glm::vec3 edge2 = glm::vec3(nextVert->pos - curVertex->pos);

            glm::vec4 normal = glm::vec4(glm::normalize(glm::cross(edge1, edge2)), 0.0f);

            if(fequal(glm::length(normal), 0.0f))
            {
                normal = faces[i]->normal;
                continue;
            }

            nor.push_back(normal);
            pos.push_back(curVertex->pos);
            //int n = 0;
            //col.push_back(glm::vec4(curVertex->weightMap[n], curVertex->weightMap[n], curVertex->weightMap[n], 1.0f));
            col.push_back(faces[i]->color);
            first = first->nextEdge;
        }while(first != faces[i] -> edge);
        unsigned int postPosCount = pos.size();
        for(int j = 0; j < postPosCount - prePosCount - 2; j++)
        {
            idx.push_back(prePosCount);
            idx.push_back(prePosCount + j + 1);
            idx.push_back(prePosCount + j + 2);
            count += 3;
        }
    }
}

void Mesh::update()
{
    // Set data vector
    std::vector<glm::vec4> pos;
    std::vector<glm::vec4> col;
    std::vector<glm::vec4> nor;
    std::vector<unsigned int> idx;
    std::vector<glm::ivec4> jointIndex; // Joint indices
    std::vector<glm::vec4> jointMaxWeight; // Joint weights

    count = 0;

    for(int i = 0; i <  faces.size(); i++)
    {
        HalfEdge* first = faces[i]->edge;
        unsigned int prePosCount = pos.size();
        do {
            Vertex* curVertex = first->nextVertex;

            std::vector<float> jointWeight;
            std::vector<float> jointSumWeight;
            std::vector<int> jointWeightIndices;
            std::unordered_map<int, bool> maxMap;

            Vertex* lastVert = first->symEdge->nextVertex;
            Vertex* nextVert = first->nextEdge->nextVertex;
            glm::vec3 edge1 = glm::vec3(curVertex->pos - lastVert->pos);
            glm::vec3 edge2 = glm::vec3(nextVert->pos - curVertex->pos);

            glm::vec4 normal = glm::vec4(glm::normalize(glm::cross(edge1, edge2)), 0.0f);

            if(std::isnan(normal.x))
            {
                normal = faces[i]->normal;
            }

            nor.push_back(normal);
            pos.push_back(curVertex->pos);
            //int n = 1;
            //float color = curVertex->weightMap[n];
            //col.push_back(glm::vec4(curVertex->weightMap[0], curVertex->weightMap[1], curVertex->weightMap[13], 1.0f));
            col.push_back(faces[i]->color);

            first = first->nextEdge;
        }while(first != faces[i] -> edge);

        unsigned int postPosCount = pos.size();
        for(int j = 0; j < postPosCount - prePosCount - 2; j++)
        {
            idx.push_back(prePosCount);
            idx.push_back(prePosCount + j + 1);
            idx.push_back(prePosCount + j + 2);
            count += 3;
        }
    }
}

// Construct a mesh from file
void Mesh::constructFileMesh(std::string& fileName)
{
    srand(time(NULL));
    // Clear all mesh components
    vertices.clear();
    halfEdges.clear();
    faces.clear();

    faceVerticesCount = 0;

    // Index vector
    std::vector<int> indices;
    std::unordered_map<std::string, HalfEdge*> symMap;

    // Open obj file
    const char* cFileName = (char*)fileName.c_str();
    std::ifstream fileIn;
    fileIn.open(cFileName);

    // Read file
    std::string lineBuffer = "";
    while(std::getline(fileIn, lineBuffer))
    {
        std::istringstream iss(lineBuffer);
        std::string split;

        std::vector<std::string> lineSegment;

        bool first = true;

        if(lineBuffer == "")
        {
            continue;
        }

        // Split line segment
        while(std::getline(iss, split, ' '))
        {
           lineSegment.push_back(split);
           if(first)
           {
               first = false;
               if(split != "v" && split != "f")
               {
                   break;
               }
           }
        }

        if(lineSegment[0] == "v")
        {
            // Read vertex info
            float x = float(atof((char*)lineSegment[1].c_str()));
            float y = float(atof((char*)lineSegment[2].c_str()));
            float z = float(atof((char*)lineSegment[3].c_str()));

            std::unique_ptr<Vertex> curVertex = std::make_unique<Vertex>();
            glm::vec4 pos(x, y, z, 1.0f);
            curVertex->pos = pos;

            vertices.push_back(std::move(curVertex));
        }
        else if(lineSegment[0] == "f")
        {
            faceVerticesCount = lineSegment.size() - 1;
            // Read vertex index
            for(int i = 1; i < lineSegment.size(); i++)
            {
                std::string indexString;
                std::istringstream indexReader(lineSegment[i]);
                if(std::getline(indexReader, indexString, '/'))
                {
                    indices.push_back(atoi(indexString.c_str()));
                }

            }
        }
    }

    // Close file
    fileIn.close();

    // Set mesh data structure
    for(int i = 0; i < indices.size(); i = i + faceVerticesCount)
    {
        // Set a new face
        std::unique_ptr<Face> f = std::make_unique<Face>();

        // Get vertices
        std::vector<Vertex*> faceVertices;
        for(int j = 0; j < faceVerticesCount; j++)
        {
            faceVertices.push_back(vertices[indices[i + j] - 1].get());
        }

        // Create halfEdges
        std::vector<HalfEdge*> faceHalfEdges;
        for(int k = 0; k < faceVerticesCount; k++)
        {
            std::unique_ptr<HalfEdge> he = std::make_unique<HalfEdge>();
            he->nextVertex = faceVertices[k];
            faceVertices[k]->edge = he.get();
            he->face = f.get();
            f->edge = he.get();
            faceHalfEdges.push_back(he.get());

            // Add halfEdges
            halfEdges.push_back(std::move(he));
        }

        std::string middle = "t";
        // Setup halfEdge pointers
        for(int w = 0; w < faceVerticesCount; w++)
        {
            // Set nextEdge
            if(w == faceVerticesCount - 1)
            {
                faceHalfEdges[w]->nextEdge = faceHalfEdges[0];
            }
            else
            {
                faceHalfEdges[w]->nextEdge = faceHalfEdges[w + 1];
            }
            // Set symEdge
            std::string heSymString = "";
            if(w == 0)
            {
                heSymString = std::to_string(indices[i]) + middle + std::to_string(indices[i + faceVerticesCount - 1]);
            }
            else
            {
                heSymString = std::to_string(indices[i + w]) + middle + std::to_string(indices[i + w - 1]);
            }
            // Get sym edge
            HalfEdge* heSym = symMap[heSymString];

            // Set sym edge
            if(heSym == NULL)
            {
                std::string name = "";
                if(w == 0)
                {
                    name = std::to_string(indices[i + faceVerticesCount - 1]) + middle + std::to_string(indices[i]);
                }
                else
                {
                    name = std::to_string(indices[i + w - 1]) + middle + std::to_string(indices[i + w]);
                }
                symMap[name] = faceHalfEdges[w];
            }
            else
            {
                faceHalfEdges[w]->symEdge = heSym;
                heSym->symEdge = faceHalfEdges[w];
                symMap.erase(heSymString);
            }
        }

        // Set face color

        float r = rand() / float(RAND_MAX);
        float g = rand() / float(RAND_MAX);
        float b = rand() / float(RAND_MAX);

        f->color = glm::vec4(r, g, b, 1.0f);
        faces.push_back(std::move(f));
    }


    for (int i = 0; i < vertices.size(); i++)
    {
        vertices.at(i)->edge = vertices.at(i)->edge->symEdge;
    }
}

Mesh::~Mesh() {}


