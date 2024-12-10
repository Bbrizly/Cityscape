#pragma once
#include <vector>
#include <cmath>
#include "GeometryUtils.h"
#include <random>
#include "Vertex.h"
#include "GeometryUtils.h"
#include <glm/glm.hpp>

class PolygonUtils {
public:
    static double distanceBetweenPoints(const Point& p1, const Point& p2);
    static std::pair<Point,Point> findLargestEdge(std::vector<Point> polygon);
    static std::pair<std::vector<Point>, std::vector<Point>> clipPolygon(std::vector<Point> polygon, Line l);
    static Point getCentroid(std::vector<Point> polygon);
    static bool isConvex(const std::vector<Point>& polygon);
    static Line moveLineToCenter(Line l, std::vector<Point> polygon);
    static std::vector<Point> scalePolygon(const std::vector<Point>& polygon, float scaleFactor);
    static float calculatePolygonArea(const std::vector<Point>& polygon);
    static bool findSmallestEdgeAmount(std::vector<Point> polygon, float minEdgeLength);
    static std::vector<Vertex> fanTriangulatePolygon(const std::vector<Point>& polygon, 
                        const glm::vec3& normal, float height, float layer,
                        GLubyte r, GLubyte g, GLubyte b, GLubyte a);
};
