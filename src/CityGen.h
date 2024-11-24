#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Cube.h"
#include "Vertex.h"
#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
using namespace std;

struct Point {double x, y;};

struct Line {double a, b, c;}; // Line equation: ax + by + c = 0

struct Edge {Point start, end;};

class CityGen {
private:

    GLuint m_program;
    GLuint m_vao, m_vbo;
    vector<Vertex> m_vertices;
    int m_numVertices;

    // Voronoi parameters
    const int numSites = 11;
    const double minX = -300.0;
    const double maxX = 300.0;
    const double minY = -200.0;
    const double maxY = 200.0;
    const double epsilon = 1e-9;

    vector<Point> m_sites;
    vector<vector<Point>> m_voronoiCells;
    vector<vector<Point>> m_chunks;

    // Functions for Voronoi diagram
    vector<Point> findIntersectionsWithBoundary(const Line& line);
    void CreateBuildingsAlongLine(const Point& start, const Point& end, int numBuildings);
    vector<Point> generateSites(int numSites, unsigned int seed);
    Line perpendicularBisector(const Point& p1, const Point& p2);
    double evaluate(const Line& l, const Point& p);
    bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection);
    vector<Point> clipPolygon(const vector<Point>& polygon, const Line& l, bool keepPositiveSide);
    void computeVoronoiDiagram();
    void computeChunks();
    void buildVertexData();
    
public:
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const glm::mat4& proj, const glm::mat4& view, float m_timer);
};