#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Cube.h"
#include "Vertex.h"
#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
using namespace std;
using namespace glm;

struct Point {double x, y;};

struct Line {double a, b, c;}; // Line equation: ax + by + c = 0

struct Edge {Point start, end;};

class CityGen {
private:

    vector<vector<Point>> debugs;
    vector<vector<Point>> debugs1;
    vector<vector<Point>> debugs2;
    GLuint otherShader;

    GLuint m_program;
    GLuint m_vao, m_vbo;
    vector<Vertex> m_vertices;
    int m_numVertices;

    // Voronoi parameters
    const int numSites = 10;
    const double minX = -200.0;
    const double maxX = 200.0;
    const double minY = -100.0;
    const double maxY = 100.0;
    const double epsilon = 1e-9;

    vector<Point> m_sites;
    vector<vector<Point>> m_voronoiCells;
    vector<vector<Point>> m_chunks;
    vector<vector<Point>> m_blocks;
    vector<vector<Point>> m_buildings;
    Point getDirectionVector(const Point& from, const Point& to);
    
    pair<pair<Point, Point>, vec2> getEdgeWithInwardDirection(
    vector<Point>& polygon,
    size_t edgeIndex);
    
    Point getCentroid(vector<Point> polygon);

    Line moveLineInDirection(const Line& line, const Point& direction, double distance);

    bool findSmallestEdgeAmount(vector<Point> polygon, float minEdgeLength);

    bool isPolygonClockwise(const vector<Point>& polygon);
    pair<Point,Point> findLargestEdge(vector<Point> polygon);
    // Line findLargestEdge(vector<Point> polygon);
    float findSmallestEdgeAmount(vector<Point> polygon);
    Point findMidpoint(const Point& p1, const Point& p2);
    void CreateCube(Point x, float height);
    Line moveToPoint(Line l, Point p);
    Line moveLineToCenter(Line l, vector<Point> polygon);

    void addDebug1(Point p);
    void addLineDebug1(Point p1,Point p2);

    void addDebug(Point p);
    void addLineDebug(Point p1,Point p2);
    Point movePointInDirection(const Point& point, const Point& direction, double distance);

    const float maximumChunkSize = 5000.0f;

    float calculatePolygonArea(const vector<Point>& polygon);

    vector<Line> generateSweepLines(const Point& origin, const Point& direction, double minSpacing, double maxSpacing, double maxDistance);

    bool isConvex(const vector<Point>& polygon);
    Point normalizeVector(double x, double y);
    Point perpendicularVector(double x, double y);
    vector<Point> offsetPolygonInward(const vector<Point>& polygon, double offsetDistance);

    vector<Point> scalePolygon(const vector<Point>& polygon, float scaleFactor);

    // Functions for Voronoi diagram
    vector<Point> findIntersectionsWithBoundary(const Line& line);

    Line CreateLineFromPoints(const Point& p1, const Point& p2);

    void CreateCubesAlongLine(const Point& start, const Point& end, int numBuildings);

    vector<Point> generateSites(int numSites, unsigned int seed);
    Line perpendicularBisector(const Point& p1, const Point& p2);
    double evaluate(const Line& l, const Point& p);
    bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection);
    // vector<Point> clipPolygon(const vector<Point>& polygon, const Line& l, bool keepPositiveSide);
    
    pair<vector<Point>, vector<Point>> clipPolygon(vector<Point> polygon, Line l);
    vector<vector<Point>> splitToBlocks(vector<Point>& polygon, vector<Line>& sweepLines);
    
    void computeVoronoiDiagram();
    void computeChunks();
    void sweepToBlocks();
    void buildVertexData();
    void buildBuildings();

    
public:
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const mat4& proj, const mat4& view, float m_timer);
};