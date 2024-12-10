#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include "Types.h"
#include <glm/glm.hpp>

class GeometryUtils {
public:
    static bool findLineIntersection(const Line& l1, const Line& l2, Point& intersection);
    static double evaluate(const Line& l, const Point& p);
    static Line CreateLineFromPoints(const Point& p1, const Point& p2);
    static Point findMidpoint(const Point& p1, const Point& p2);
    static Line perpendicularBisector(const Point& p1, const Point& p2);
    static bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l1, Point& intersection);
    static Point getDirectionVector(const Point& from, const Point& to);
    static Point getPerpendicularDirVector(std::pair<Point,Point> edge, const Point& centroid);
    static Point movePointInDirection(const Point& point, const Point& direction, double distance);
    static Line moveToPoint(Line l, Point p);
    static Point normalizeVector(double x, double y);
    static Point perpendicularVector(double x, double y);
    static Point oppositeVector(double x, double y);
    static Line makePerpendicularLine(const Line& originalLine);
    static glm::vec3 calculateQuadNormal(const Point& p1, const Point& p2);
};
