#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Vertex.h"
#include "Types.h"
#include "GeometryUtils.h"
#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include "BuildingGen.h"
#include "PolygonUtils.h"
#include "Debug.h"
#include "Voronoi.h"
#include "DrawRoad.h"
#include <functional>
#include <random>
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <stb_image.h>
using namespace std; using namespace glm;

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
    vector<vector<Point>> debugs;
    vector<vector<Point>> debugs1;
    vector<vector<Point>> debugs2;
    
    BuildingGen* m_buildingGen;
    void initBuildingGenerator();

    bool Debug = true;
    int m_numVertices;
    int m_numLines;

    //Textures layers
    float wall = 0;
    float wall2 = 1;
    float wall3 = 2;
    float windows = 3;
    float windows2 = 4;
    float windows3 = 5;
    float roof = 6;
    float sidewalk = 7;
    float asphalt = 8;

    // Voronoi diagram parameters
    const int numSites = 25;
    float moveAmount = 25.0f;
    float maxMoveAmount = 50.0f;
    float minMoveAmount = 40.0f;

    float minPolygonArea = (moveAmount * moveAmount )/4;

    float minEdge = 2.0f;
    float minEdge2 = 0.5f;

    const double minX = -500.0;
    const double maxX = 500.0;
    const double minY = -400.0;
    const double maxY = 400.0;
    float cyceLength = 48.0f;
    const double epsilon = 1e-9;

    vector<Point> m_sites;
    vector<vector<Point>> m_voronoiCells;
    vector<vector<Point>> m_chunks;
    vector<vector<Point>> m_blocks;
    vector<vector<Point>> m_strips;
    vector<Building> m_buildings;
    Point districtCenter;
    float districtRadius = 0.0f;
    
    const int ResidentialMinStories = 2;
    const int ResidentialMaxStories = 8;

    const int CommercialMinStories = 8;
    const int CommercialMaxStories = 15;

    const int IndustrialMinStories = 20;
    const int IndustrialMaxStories = 35;

    const float baseAddition0 = 0.8f;    // added to the top to make it seamlesss
    const float baseAddition2 = 0.6f;
    const float baseStoryHeight = 10.0f;    // Height per story
    //normally is 10 for 618 x 600 but im 10.8 for extra seamlessness in texturing 

    // Special Industrial Building Controls
    const float specialBuildingChance = 0.1f;           // 10% chance
    const int extraBuildingCountMin = 3;                // Minimum extra buildings
    const int extraBuildingCountMax = 4;                // Maximum extra buildings
    const float extraBuildingScaleFactor = 0.8f;         // Scale factor for extra buildings
    const float extraBuildingHeightMultiplier = 0.5f;    // Height multiplier for extra buildings

    //Data
    vector<Vertex> m_vertices;
    vector<Vertex> m_lines;
    
    void pushVertexData(wolf::VertexBuffer*& vBuffer, wolf::VertexDeclaration*& vDecl, vector<Vertex>& vertices);
    vec3 calculateQuadNormal(const Point& p1, const Point& p2);

    Building determineBuildingDetails(const vector<Point>& polygon);
    void BuildingToVerticies(const Building& building, vector<Vertex>& m_vertices, float ground);

    void buildVertexData();
    void computeChunks();
    void sweepToBlocks();

public:
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const mat4& proj, const mat4& view, float m_timer);
};
