#pragma once
#include <vector>
#include "GeometryUtils.h"

class Debug {
public:
    static void addDebug1(std::vector<std::vector<Point>>& debugs1, Point p);
    static void addLineDebug1(std::vector<std::vector<Point>>& debugs1, Point p1,Point p2);
    static void addLineToVector1(std::vector<std::vector<Point>>& debugs1, Point p1,Point vector);

    static void addDebug(std::vector<std::vector<Point>>& debugs2, Point p);
    static void addLineDebug(std::vector<std::vector<Point>>& debugs2, Point p1,Point p2);
    static void addLineToVector(std::vector<std::vector<Point>>& debugs2, Point p1,Point vector);

    static std::vector<Point> findIntersectionsWithBoundary(const Line& line, double minX, double maxX, double minY, double maxY);
};
