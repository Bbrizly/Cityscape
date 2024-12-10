#include "Debug.h"

void Debug::addDebug1(std::vector<std::vector<Point>>& debugs1, Point p)
{
    debugs1.push_back({{p.x-5.0f,p.y -5.0f}, {p.x+5.0f,p.y +5.0f}});
    debugs1.push_back({{p.x-5.0f,p.y +5.0f}, {p.x +5.0f,p.y -5.0f}});
}

void Debug::addLineDebug1(std::vector<std::vector<Point>>& debugs1, Point p1,Point p2)
{
    debugs1.push_back({p1,p2});
}

void Debug::addLineToVector1(std::vector<std::vector<Point>>& debugs1, Point p1,Point vector)
{
    addDebug1(debugs1,p1);
    addLineDebug1(debugs1,p1,{p1.x + vector.x * 10.0f, p1.y + vector.y * 10.0f});
}

void Debug::addDebug(std::vector<std::vector<Point>>& debugs2, Point p)
{
    debugs2.push_back({{p.x-5.0f,p.y -5.0f}, {p.x+5.0f,p.y +5.0f}});
    debugs2.push_back({{p.x-5.0f,p.y +5.0f}, {p.x +5.0f,p.y -5.0f}});
}

void Debug::addLineDebug(std::vector<std::vector<Point>>& debugs2, Point p1,Point p2)
{
    debugs2.push_back({p1,p2});
}

void Debug::addLineToVector(std::vector<std::vector<Point>>& debugs2, Point p1,Point vector)
{
    addDebug(debugs2,p1);
    addLineDebug(debugs2,p1,{p1.x + vector.x * 10.0f, p1.y + vector.y * 10.0f});
}

std::vector<Point> Debug::findIntersectionsWithBoundary(const Line& line, double minX, double maxX, double minY, double maxY) {
    std::vector<Point> intersections;
    Line top = {0, 1, -maxY};
    Line bottom = {0, 1, -minY};
    Line left = {1, 0, -minX};
    Line right = {1, 0, -maxX};
    Point intersection;
    if (GeometryUtils::findLineIntersection(line, top, intersection)) {
        if (intersection.x >= minX && intersection.x <= maxX) {
            intersections.push_back(intersection);
        }
    }
    if (GeometryUtils::findLineIntersection(line, bottom, intersection)) {
        if (intersection.x >= minX && intersection.x <= maxX) {
            intersections.push_back(intersection);
        }
    }
    if (GeometryUtils::findLineIntersection(line, left, intersection)) {
        if (intersection.y >= minY && intersection.y <= maxY) {
            intersections.push_back(intersection);
        }
    }
    if (GeometryUtils::findLineIntersection(line, right, intersection)) {
        if (intersection.y >= minY && intersection.y <= maxY) {
            intersections.push_back(intersection);
        }
    }
    return intersections;
}
