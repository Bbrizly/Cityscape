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

// Basic geometric structures
struct Point { double x, y; };
struct Line { double a, b, c; }; // Line equation: ax + by + c = 0
struct Edge { Point start, end; };

// CityGen Class Definition
class CityGen {
private:
    // Debugging vectors
    vector<vector<Point>> debugs;
    vector<vector<Point>> debugs1;
    vector<vector<Point>> debugs2;

    // OpenGL related variables
    GLuint otherShader;
    GLuint m_program;
    GLuint textureID;
    GLuint m_vao, m_vbo;
    GLuint m_vaoLines, m_vboLines;
    vector<Vertex> m_vertices;
    vector<Vertex> m_lines;
    bool Debug = false;
    int m_numVertices;
    int m_numLines;

    // Voronoi diagram parameters
    const int numSites = 23;
    float moveAmount = 25.0f;
    float maxMoveAmount = 30.0f;
    float minMoveAmount = 20.0f;
    float minPolygonArea = (moveAmount * moveAmount )/4;//50.0f;
    float minEdge = 2.0f;
    float minEdge2 = 0.5f;

    const double minX = -500.0;
    const double maxX = 500.0;
    const double minY = -400.0;
    const double maxY = 400.0;
    const double epsilon = 1e-9;

    // Geometric and utility functions
    GLuint loadTexture(const std::string& filepath);
    void drawCrosswalk(Point base, Point direction, float lineLength, float spacing, float lineHeight);
    vector<pair<Point, Point>> findSharedEdges(const vector<vector<Point>>& polygons);
    void addRoadDecals(Point p1, Point p2);
    double distanceBetweenPoints(const Point& p1, const Point& p2);
    Line CreateLineFromPoints(const Point& p1, const Point& p2);
    Point findMidpoint(const Point& p1, const Point& p2);
    Line perpendicularBisector(const Point& p1, const Point& p2);
    bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection);
    double evaluate(const Line& l, const Point& p);
    Point getPerpendicularDirVector(pair<Point,Point> edge, const Point& centroid);
    Point getDirectionVector(const Point& from, const Point& to);
    Line moveLineInDirection(const Line& line, const Point& direction, double distance);
    Point movePointInDirection(const Point& point, const Point& direction, double distance);
    Line moveToPoint(Line l, Point p);
    Point normalizeVector(double x, double y);
    Point perpendicularVector(double x, double y);
    Point oppositeVector(double x, double y);
    Line makePerpendicularLine(const Line& originalLine);
    vector<Point> findIntersectionsWithBoundary(const Line& line);
    pair<Point, Point> findLargestEdge(vector<Point> polygon);
    pair<vector<Point>, vector<Point>> clipPolygon(vector<Point> polygon, Line l);
    Point getCentroid(vector<Point> polygon);
    bool isConvex(const vector<Point>& polygon);
    Line moveLineToCenter(Line l, vector<Point> polygon);
    vector<Point> scalePolygon(const vector<Point>& polygon, float scaleFactor);
    float calculatePolygonArea(const vector<Point>& polygon);
    bool findSmallestEdgeAmount(vector<Point> polygon, float minEdgeLength);

    // Voronoi and chunk processing
    vector<Point> generateSites(int numSites, unsigned int seed);
    void computeVoronoiDiagram();
    void computeChunks();
    void sweepToBlocks();
    void buildVertexData();

    // Debugging functions
    void addDebug1(Point p);
    void addLineDebug1(Point p1, Point p2);
    void addDebug(Point p);
    void addLineDebug(Point p1, Point p2);
    void addLineToVector(Point p1, Point vector);
    void addLineToVector1(Point p1, Point vector);

    // Data structures for Voronoi and building generation
    vector<Point> m_sites;
    vector<vector<Point>> m_voronoiCells;
    vector<vector<Point>> m_chunks;
    vector<vector<Point>> m_blocks;
    vector<vector<Point>> m_strips;
    vector<vector<Point>> m_buildings;

    // Constants
    const float maximumChunkSize = 5000.0f;

public:
    // Public interface
    void generate(GLuint program);
    void deGenerate();
    void reGenerate();
    void render(const mat4& proj, const mat4& view, float m_timer);
};
