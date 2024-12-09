#pragma once
#include "GeometryUtils.h"
#include "Vertex.h"
#include <functional>
#include <vector>

class DrawRoad {
public:
    static void drawCrosswalk(std::vector<Vertex>& m_lines, Point base, Point direction, float lineLength, float spacing, float lineHeight, float sidewalk);
    static void addRoadDecals(std::vector<Vertex>& m_lines, Point p1, Point p2, 
                              float sidewalk, 
                              std::function<double(const Point&, const Point&)> distanceBetweenPoints, 
                              std::function<Point(const Point&, const Point&)> getDirectionVector,
                              std::function<Point(const Point&, const Point&, double)> movePointInDirection);
};
