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
    const int numSites = 2;
    const double minX = -100.0;
    const double maxX = 100.0;
    const double minY = -100.0;
    const double maxY = 100.0;
    // const double minX = 0;
    // const double maxX = 200;
    // const double minY = 0;
    // const double maxY = 200;
    const double epsilon = 1e-9;

    vector<Point> m_sites;
    vector<vector<Point>> m_voronoiCells;

    // Functions for Voronoi diagram
    vector<Point> findIntersectionsWithBoundary(const Line& line);
    void CreateBuildingsAlongLine(const Point& start, const Point& end, int numBuildings);
    vector<Point> generateSites(int numSites, unsigned int seed);
    Line perpendicularBisector(const Point& p1, const Point& p2);
    double evaluate(const Line& l, const Point& p);
    bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection);
    vector<Point> clipPolygon(const vector<Point>& polygon, const Line& l);
    void computeVoronoiDiagram();
    void buildVertexData();
    
public:
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const glm::mat4& proj, const glm::mat4& view, float m_timer);
};

/*
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
    const int numSites = 15;
    const double minX = -100.0;
    const double maxX = 100.0;
    const double minY = -100.0;
    const double maxY = 100.0;
    const double epsilon = 1e-9;

    vector<Point> m_sites;
    vector<vector<Point>> m_voronoiCells;

    // Functions for Voronoi diagram
    void CreateBuildingsAlongLine(const Point& start, const Point& end, int numBuildings);
    vector<Point> generateSites(int numSites, unsigned int seed);
    Line perpendicularBisector(const Point& p1, const Point& p2);
    double evaluate(const Line& l, const Point& p);
    bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection);
    vector<Point> clipPolygon(const vector<Point>& polygon, const Line& l, const Point& site);
    void computeVoronoiDiagram();
    void buildVertexData();
    
public:
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const glm::mat4& proj, const glm::mat4& view, float m_timer);
};

*/
/*
class CityGen {
private:
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_program = 0;

public:

    // struct Site {float x, z;};

    // void computeCatmullRomSpline(const std::vector<Site>& sites, int pointsPerSegment);
    // void generateRoadMesh(const std::vector<glm::vec3>& splinePoints, float roadWidth);
    

    void generate(GLuint program);
    void reGenerate();
    void deGenerate();
    void render(const glm::mat4& proj, const glm::mat4& view, float m_timer);
};
*/