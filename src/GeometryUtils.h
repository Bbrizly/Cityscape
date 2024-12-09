#pragma once
#include <vector>
#include <cmath>
#include <iostream>

struct Point { double x, y; };
struct Line { double a, b, c; };

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
    static Line moveLineInDirection(const Line& line, const Point& direction, double distance);
    static Point movePointInDirection(const Point& point, const Point& direction, double distance);
    static Line moveToPoint(Line l, Point p);
    static Point normalizeVector(double x, double y);
    static Point perpendicularVector(double x, double y);
    static Point oppositeVector(double x, double y);
    static Line makePerpendicularLine(const Line& originalLine);
};
