#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Vertex.h"
#include "GeometryUtils.h"
#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>

// Forward declarations of utility classes
class PolygonUtils;
class Voronoi;
class Debug;
class DrawRoad;

class CityGen {
private:
    wolf::Program* m_program;
    wolf::Texture* m_buildingTexture;
    wolf::Texture* m_arrayTexture;
    
    wolf::VertexBuffer* m_vertexBuffer;
    wolf::VertexBuffer* m_lineBuffer;
    wolf::VertexDeclaration* m_vertexDecl;
    wolf::VertexDeclaration* m_lineDecl;

    // Debugging vectors
    std::vector<std::vector<Point>> debugs;
    std::vector<std::vector<Point>> debugs1;
    std::vector<std::vector<Point>> debugs2;

    bool Debug = true;
    int m_numVertices;
    int m_numLines;

    //Textures layers
    float wall = 0;
    float wall2 = 1;
    float roof = 2;
    float sidewalk = 3;
    float windows = 4;
    float windows2 = 5;
    float asphalt = 6;

    // Voronoi diagram parameters
    const int numSites = 125;
    float moveAmount = 25.0f;
    float maxMoveAmount = 30.0f;
    float minMoveAmount = 20.0f;
    float minPolygonArea = (moveAmount * moveAmount )/4;
    float minEdge = 2.0f;
    float minEdge2 = 0.5f;

    const double minX = -1500.0;
    const double maxX = 1500.0;
    const double minY = -1400.0;
    const double maxY = 1400.0;
    const double epsilon = 1e-9;

    std::vector<Point> m_sites;
    std::vector<std::vector<Point>> m_voronoiCells;
    std::vector<std::vector<Point>> m_chunks;
    std::vector<std::vector<Point>> m_blocks;
    std::vector<std::vector<Point>> m_strips;
    std::vector<std::vector<Point>> m_buildings;

    //Data
    std::vector<Vertex> m_vertices;
    std::vector<Vertex> m_lines;

    wolf::Texture* CreateArrayTextureFromFiles(const std::vector<std::string>& filePaths);
    void pushVertexData(wolf::VertexBuffer*& vBuffer, wolf::VertexDeclaration*& vDecl, std::vector<Vertex>& vertices);
    glm::vec3 calculateQuadNormal(const Point& p1, const Point& p2);

    void buildVertexData();
    void computeChunks();
    void sweepToBlocks();

public:
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const glm::mat4& proj, const glm::mat4& view, float m_timer);
};
