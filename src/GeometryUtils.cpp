#include "GeometryUtils.h"

bool GeometryUtils::findLineIntersection(const Line& l1, const Line& l2, Point& intersection) {
    double det = l1.a * l2.b - l2.a * l1.b;
    if (std::abs(det) < 1e-9) {
        return false;
    }
    intersection.x = (l2.b * -l1.c - l1.b * -l2.c) / det;
    intersection.y = (l1.a * -l2.c - l2.a * -l1.c) / det;
    return true;
}

double GeometryUtils::evaluate(const Line& l, const Point& p) {
    return (l.a * p.x + l.b * p.y + l.c);
}

Line GeometryUtils::CreateLineFromPoints(const Point& p1, const Point& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    if (std::abs(dx) < 1e-9) {
        return {1.0, 0.0, -p1.x};
    }
    double slope = dy / dx;
    double a = -slope;
    double b = 1.0;
    double c = -(a * p1.x + b * p2.y);
    return {a, b, c};
}

Point GeometryUtils::findMidpoint(const Point& p1, const Point& p2) {
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;
    return {mx,my};
}

Line GeometryUtils::perpendicularBisector(const Point& p1, const Point& p2) {
    Point midpoint = findMidpoint(p1,p2);
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    if (std::abs(dx) < 1e-9) {
        return {0.0, 1.0, -(midpoint.y)};
    }
    double slope = dy / dx;
    double perpSlope = -dx / dy;
    double a = -perpSlope;
    double b = 1.0;
    double c = -(a * midpoint.x + b * midpoint.y);
    if ((a * p1.x + b * p1.y + c) < 0) {
        a = -a; b = -b; c = -c;
    }
    return {a, b, c};
}

bool GeometryUtils::lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l1, Point& intersection) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    if (std::abs(dx) < 1e-9) {
        Line l2 = {1.0, 0.0, -p1.x};
        bool ret = findLineIntersection(l1, l2, intersection);
        if (ret) {
            double minY = std::min(p1.y, p2.y) - 1e-9;
            double maxY = std::max(p1.y, p2.y) + 1e-9;
            if (intersection.y >= minY && intersection.y <= maxY) {
                return true;
            }
        }
        return false;
    }
    double slope = dy / dx;
    Line l2 = {-slope, 1.0, slope * p1.x - p1.y};
    bool ret = findLineIntersection(l1, l2, intersection);

    if (ret) {
        double minX = std::min(p1.x, p2.x) - 1e-9;
        double maxX = std::max(p1.x, p2.x) + 1e-9;
        double minY = std::min(p1.y, p2.y) - 1e-9;
        double maxY = std::max(p1.y, p2.y) + 1e-9;
        if (intersection.x >= minX && intersection.x <= maxX &&
            intersection.y >= minY && intersection.y <= maxY) {
            return true;
        }
    }
    return false;
}

Point GeometryUtils::getDirectionVector(const Point& from, const Point& to) {
    double dx = to.x - from.x;
    double dy = to.y - from.y;
    double magnitude = sqrt(dx * dx + dy * dy);
    if (magnitude < 1e-9) {
        std::cerr << "Error: Zero-length direction vector." << std::endl;
        return {0.0, 0.0};
    }
    return {dx / magnitude, dy / magnitude};
}

Point GeometryUtils::getPerpendicularDirVector(std::pair<Point,Point> edge, const Point& centroid) {
    Point edgeMidpoint = findMidpoint(edge.first,edge.second);
    Point edgeDirectionVector = getDirectionVector(edgeMidpoint,centroid);
    double dx = edge.second.x - edge.first.x;
    double dy = edge.second.y - edge.first.y;
    Point perp1 = perpendicularVector(dx, dy);
    Point perp2 = oppositeVector(perp1.x, perp1.y);
    if ((perp1.x * edgeDirectionVector.x + perp1.y * edgeDirectionVector.y) > 0) {
        return normalizeVector(perp1.x, perp1.y);
    } else {
        return normalizeVector(perp2.x, perp2.y);
    }
}

Line GeometryUtils::moveLineInDirection(const Line& line, const Point& direction, double distance) {
    double displacement = line.a * direction.x * distance + line.b * direction.y * distance;
    double newC = line.c + displacement;
    return {line.a, line.b, newC};
}

Point GeometryUtils::movePointInDirection(const Point& point, const Point& direction, double distance) {
    double magnitude = sqrt(direction.x * direction.x + direction.y * direction.y);
    if (magnitude < 1e-9) {
        std::cerr << "Error: Zero-length direction vector." << std::endl;
        return point;
    }
    Point normalizedDirection = {direction.x / magnitude, direction.y / magnitude};
    double newX = point.x + normalizedDirection.x * distance;
    double newY = point.y + normalizedDirection.y * distance;
    return {newX, newY};
}

Line GeometryUtils::moveToPoint(Line l, Point p) {
    l.c = -(l.a * p.x + l.b * p.y);
    return l;
}

Point GeometryUtils::normalizeVector(double x, double y) {
    double mag = std::sqrt(x * x + y * y);
    if (mag == 0.0) return {0.0, 0.0};
    return {x / mag, y / mag};
}

Point GeometryUtils::perpendicularVector(double x, double y) {
    return {-y, x};
}

Point GeometryUtils::oppositeVector(double x, double y) {
    return {-x, -y};
}

Line GeometryUtils::makePerpendicularLine(const Line& originalLine) {
    double a = originalLine.b;
    double b = -originalLine.a;
    double c = originalLine.c;
    return {a, b, c};
}
