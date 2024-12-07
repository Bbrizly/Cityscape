#include "CityGen.h"
#include "Vertex.h"
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <iostream>
#include <glm/glm.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stb_image.h>
using namespace std;
using namespace glm;
#pragma region ALL
#pragma region Mathematical shit
bool findLineIntersection(const Line& l1, const Line& l2, Point& intersection) {
    // Calculate the determinant
    double det = l1.a * l2.b - l2.a * l1.b;

    if (std::abs(det) < 1e-9) {
        // Lines are parallel or coincident
        return false;
    }

    // Use Cramer's rule to solve for the intersection
    intersection.x = (l2.b * -l1.c - l1.b * -l2.c) / det;
    intersection.y = (l1.a * -l2.c - l2.a * -l1.c) / det;
    return true;
}
double CityGen::evaluate(const Line& l, const Point& p) {
    return (l.a * p.x + l.b * p.y + l.c);
}
Line CityGen::CreateLineFromPoints(const Point& p1, const Point& p2)
{
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    if (std::abs(dx) < 1e-9) {
        // Vertical line equation: x = p1.x
        return {1.0, 0.0, -p1.x};
    }

    double slope = dy / dx;

    double a = -slope;
    double b = 1.0;
    double c = -(a * p1.x + b * p2.y);

    return {a, b, c};
}
Point CityGen::findMidpoint(const Point& p1, const Point& p2)
{
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;
    return {mx,my};
}
Line CityGen::perpendicularBisector(const Point& p1, const Point& p2) {

    Point midpoint = findMidpoint(p1,p2);

    // Slope of the line between p1 and p2
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    if (std::abs(dx) < 1e-9) {
        // Vertical line; perpendicular bisector is horizontal
        return {0.0, 1.0, -(midpoint.y)};
    }

    double slope = dy / dx;
    // Slope of the perpendicular bisector
    double perpSlope = -dx / dy;

    // Line equation: (y - my) = perpSlope * (x - mx)
    // Rearranged to a*x + b*y + c = 0
    double a = -perpSlope;
    double b = 1.0;
    double c = -(a * midpoint.x + b * midpoint.y);

    // Ensure the normal vector points towards p1
    if ((a * p1.x + b * p1.y + c) < 0) {
        a = -a;
        b = -b;
        c = -c;
    }

    return {a, b, c};
}
bool CityGen::lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l1, Point& intersection) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    // vertical lines
    if (std::abs(dx) < 1e-9) {
        // Vertical line segment; line equation x = p1.x
        Line l2;
        l2.a = 1.0;
        l2.b = 0.0;
        l2.c = -p1.x;
        bool ret = findLineIntersection(l1, l2, intersection);
        if (ret) {
            // Check if y is within the segment
            double minY = std::min(p1.y, p2.y) - 1e-9;
            double maxY = std::max(p1.y, p2.y) + 1e-9;
            if (intersection.y >= minY && intersection.y <= maxY) {
                return true;
            }
        }
        return false;
    }

    // Handle non-vertical lines
    double slope = dy / dx;
    Line l2;
    l2.a = -slope;
    l2.b = 1.0;
    l2.c = slope * p1.x - p1.y;
    bool ret = findLineIntersection(l1, l2, intersection);

    if (ret) {
        // Check if intersection lies within the segment
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
Point CityGen::getDirectionVector(const Point& from, const Point& to) {
    // Calculate the raw vector components
    double dx = to.x - from.x;
    double dy = to.y - from.y;

    // Compute the magnitude of the vector
    double magnitude = sqrt(dx * dx + dy * dy);

    // Check for zero-length vector
    if (magnitude < 1e-9) {
        cerr << "Error: Zero-length direction vector." << endl;
        return {0.0, 0.0}; // Return a zero vector in case of error
    }

    // Normalize the vector
    return {dx / magnitude, dy / magnitude};
}
Point CityGen::getPerpendicularDirVector(pair<Point,Point> edge, const Point& centroid) {
    
    Point edgeMidpoint = findMidpoint(edge.first,edge.second);
    Point edgeDirectionVector = getDirectionVector(edgeMidpoint,centroid);

    double dx = edge.second.x - edge.first.x;
    double dy = edge.second.y - edge.first.y;
    Point perp1 = perpendicularVector(dx, dy);
    Point perp2 = oppositeVector(perp1.x, perp1.y);

    // Step 4: Check which perpendicular vector points inward
    if ((perp1.x * edgeDirectionVector.x + perp1.y * edgeDirectionVector.y) > 0) {
        // Perpendicular vector is inward
        return normalizeVector(perp1.x, perp1.y);
    } else {
        // Opposite perpendicular vector is inward
        return normalizeVector(perp2.x, perp2.y);
    }
}
//     double magnitude = sqrt(line.a * line.a + line.b * line.b);
//     double newC = line.c + distance * magnitude; // Shift the line by the given distance
//     return {line.a, line.b, newC};
// }
Line CityGen::moveLineInDirection(const Line& line, const Point& direction, double distance) {
    // Adjust the line's c value by the projection of the displacement onto the normal vector (a, b)
    double displacement = line.a * direction.x * distance + line.b * direction.y * distance;
    double newC = line.c + displacement;

    // Return the moved line
    return {line.a, line.b, newC};
}
Point CityGen::movePointInDirection(const Point& point, const Point& direction, double distance) {
    // Normalize the direction vector
    double magnitude = sqrt(direction.x * direction.x + direction.y * direction.y);
    if (magnitude < 1e-9) {
        cerr << "Error: Zero-length direction vector." << endl;
        return point; // Return the original point if the direction vector is invalid
    }

    Point normalizedDirection = {direction.x / magnitude, direction.y / magnitude};

    // Compute the new position by moving in the direction
    double newX = point.x + normalizedDirection.x * distance;
    double newY = point.y + normalizedDirection.y * distance;

    return {newX, newY};
}
Line CityGen::moveToPoint(Line l, Point p) {
    l.c = -(l.a * p.x + l.b * p.y);
    return l;
}
Point CityGen::normalizeVector(double x, double y) {
    double mag = std::sqrt(x * x + y * y);
    if (mag == 0.0) return {0.0, 0.0};
    return {x / mag, y / mag};
}
Point CityGen::perpendicularVector(double x, double y) {
    return {-y, x};
}
Point CityGen::oppositeVector(double x, double y) {
    return {-x, -y};
}
Line CityGen::makePerpendicularLine(const Line& originalLine) {
    // Swap and negate the coefficients of the line
    double a = originalLine.b;
    double b = -originalLine.a;
    double c = originalLine.c; // Keep the same constant term

    return {a, b, c};
}

#pragma endregion

#pragma region Polygonn
double CityGen::distanceBetweenPoints(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
}
pair<Point,Point> CityGen::findLargestEdge(vector<Point> polygon)
{
    if (polygon.size() < 2) {
        std::cerr << "Polygon must have at least two points to form an edge." << std::endl;
        return {{0,0},{0,0}};
    }

    double maxLength = 0.0;
    std::pair<Point, Point> largestEdge = {polygon[0], polygon[0]};

    size_t n = polygon.size();
    for (size_t i = 0; i < n; ++i) {
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % n]; // Wrap around to the first point

        // Calculate the distance between p1 and p2
        double length = distanceBetweenPoints(p1,p2);

        // Update the largest edge if this one is longer
        if (length > maxLength) {
            maxLength = length;
            largestEdge = {p1, p2};
        }
    }
    
    // debugs2.push_back({largestEdge.first,largestEdge.second});
    return largestEdge;

}
// Line CityGen::findLargestEdge(vector<Point> polygon)
// {
//     if (polygon.size() < 2) {
//         std::cerr << "Polygon must have at least two points to form an edge." << std::endl;
//         return {0.0f,0.0f,0.0f};
//     }
//     double maxLength = 0.0;
//     std::pair<Point, Point> largestEdge = {polygon[0], polygon[0]};
//     size_t n = polygon.size();
//     for (size_t i = 0; i < n; ++i) {
//         const Point& p1 = polygon[i];
//         const Point& p2 = polygon[(i + 1) % n]; // Wrap around to the first point
//         // Calculate the distance between p1 and p2
//         double length = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
//         // Update the largest edge if this one is longer
//         if (length > maxLength) {
//             maxLength = length;
//             largestEdge = {p1, p2};
//         }
//     }
//     debugs2.push_back({largestEdge.first,largestEdge.second});
//     return {CreateLineFromPoints(largestEdge.first,largestEdge.second)};
// }
/*Debug pair<vector<Point>, vector<Point>> CityGen::clipPolygon(vector<Point> polygon, Line l) 
{
    vector<Point> polygonPositive;
    vector<Point> polygonNegative;

    // Debug: Function entry with input polygon and line
    cout << "Entering clipPolygon function." << endl;
    cout << "Input Polygon Points: ";
    for (const auto& p : polygon) {
        cout << "(" << p.x << ", " << p.y << ") ";
    }
    cout << endl;
    cout << "Clipping Line: " << l.a << "x + " << l.b << "y + " << l.c << " = 0" << endl;

    if (polygon.empty()) {
        cout << "Debug: The input polygon is empty. Exiting clipPolygon." << endl;
        return {polygonPositive, polygonNegative};
    }

    Point prev = polygon.back();
    double prevEval = evaluate(l, prev);

    // Debug: Initial previous point and its evaluation
    cout << "Debug: Initial previous point: (" << prev.x << ", " << prev.y << ")" << endl;
    cout << "Debug: Evaluation of previous point: " << prevEval << endl;

    for (size_t i = 0; i < polygon.size(); ++i) {
        const Point& curr = polygon[i];
        double currEval = evaluate(l, curr);

        // Determine which side each point is on
        bool prevInside = prevEval >= 0;
        bool currInside = currEval >= 0;

        // Debug: Current point and its evaluation
        cout << "Debug: Processing Point " << i + 1 << ": (" << curr.x << ", " << curr.y << ")" << endl;
        cout << "Debug: Evaluation of current point: " << currEval 
             << " (" << (currInside ? "Inside" : "Outside") << ")" << endl;

        if (currInside) {
            if (!prevInside) {
                // Edge crosses from negative to positive
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    polygonPositive.push_back(intersection);
                    polygonNegative.push_back(intersection);
                    // Debug: Intersection point when crossing from negative to positive
                    cout << "Debug: Intersection (Negative to Positive) at (" 
                         << intersection.x << ", " << intersection.y << ")" << endl;
                }
            }
            // Current point is inside positive side
            polygonPositive.push_back(curr);
            // Debug: Adding current point to positive polygon
            cout << "Debug: Added current point to positive polygon." << endl;
        }
        else {
            if (prevInside) {
                // Edge crosses from positive to negative
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    polygonPositive.push_back(intersection);
                    polygonNegative.push_back(intersection);
                    // Debug: Intersection point when crossing from positive to negative
                    cout << "Debug: Intersection (Positive to Negative) at (" 
                         << intersection.x << ", " << intersection.y << ")" << endl;
                }
            }
            // Current point is inside negative side
            polygonNegative.push_back(curr);
            // Debug: Adding current point to negative polygon
            cout << "Debug: Added current point to negative polygon." << endl;
        }

        // Update previous point and its evaluation for next iteration
        prev = curr;
        prevEval = currEval;
    }

    // Debug: Summary of clipped polygons
    cout << "Debug: Clipping completed." << endl;
    cout << "Debug: Positive Polygon has " << polygonPositive.size() << " points." << endl;
    cout << "Debug: Negative Polygon has " << polygonNegative.size() << " points." << endl;

    return {polygonPositive, polygonNegative};
}
*/
pair<vector<Point>, vector<Point>> CityGen::clipPolygon(vector<Point> polygon, Line l) 
{
    vector<Point> polygonPositive;
    vector<Point> polygonNegative;
    
    if (polygon.empty()) {
        cout << "\nPolygon is empty.";
        return {polygonPositive, polygonNegative};
    }

    Point prev = polygon.back();
    double prevEval = evaluate(l, prev);

    for (const Point& curr : polygon) {
        double currEval = evaluate(l, curr);

        // Determine which side each point is on
        bool prevInside = prevEval >= 0;
        bool currInside = currEval >= 0;

        if (currInside) {
            if (!prevInside) {
                // Edge crosses from negative to positive
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    polygonPositive.push_back(intersection);
                    polygonNegative.push_back(intersection);
                }
            }
            // Current point is inside positive side
            polygonPositive.push_back(curr);
        }
        else {
            if (prevInside) {
                // Edge crosses from positive to negative
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    polygonPositive.push_back(intersection);
                    polygonNegative.push_back(intersection);
                }
            }
            // Current point is inside negative side
            polygonNegative.push_back(curr);
        }

        prev = curr; 
        prevEval = currEval;
    }

    return {polygonPositive, polygonNegative};
}
Point CityGen::getCentroid(vector<Point> polygon)
{
    if (polygon.empty()) return {};

    // Compute centroid
    double centroidX = 0.0;
    double centroidY = 0.0;
    for (const Point& p : polygon) {
        centroidX += p.x + 0.5f;
        centroidY += p.y+ 0.5f;
    }
    centroidX /= polygon.size();
    centroidY /= polygon.size();
    return {centroidX, centroidY};
}
bool CityGen::isConvex(const std::vector<Point>& polygon) {
    if (polygon.size() < 3) return false;
    bool isConvex = true;
    bool sign = false;
    size_t n = polygon.size();

    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        size_t k = (i + 2) % n;

        double dx1 = polygon[j].x - polygon[i].x;
        double dy1 = polygon[j].y - polygon[i].y;
        double dx2 = polygon[k].x - polygon[j].x;
        double dy2 = polygon[k].y - polygon[j].y;

        double cross = dx1 * dy2 - dy1 * dx2;
        if (i == 0) {
            sign = cross > 0;
        }
        else {
            if ((cross > 0) != sign) {
                isConvex = false;
                break;
            }
        }
    }

    return isConvex;
}
Line CityGen::moveLineToCenter(Line l, vector<Point> polygon)
{
    Point centroid = getCentroid(polygon);
    return (moveToPoint(l, centroid));
}
vector<Point> CityGen::scalePolygon(const vector<Point>& polygon, float scaleFactor) {
    if (polygon.empty()) return {};

    // Compute centroid
    double centroidX = 0.0;
    double centroidY = 0.0;
    for (const Point& p : polygon) {
        centroidX += p.x;
        centroidY += p.y;
    }
    centroidX /= polygon.size();
    centroidY /= polygon.size();

    // Scale points towards centroid
    vector<Point> scaledPolygon;
    for (const Point& p : polygon) {
        double newX = centroidX + (p.x - centroidX) * scaleFactor;
        double newY = centroidY + (p.y - centroidY) * scaleFactor;
        scaledPolygon.push_back({newX, newY});
    }

    return scaledPolygon;
}
float CityGen::calculatePolygonArea(const vector<Point>& polygon) {
    float area = 0.0f;
    size_t n = polygon.size();
    for (size_t i = 0; i < n; i++) {
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % n]; // Wrap around
        area += (p1.x * p2.y - p2.x * p1.y);
    }
    return std::abs(area) / 2.0f;
}
/*float CityGen::findSmallestEdgeAmount(vector<Point> polygon)
{
    if (polygon.size() < 2) {
        cout<<"Polygon size: "<<polygon.size()<<endl;
        cerr << "SMALLESTPolygon must have at least two points to form an edge." <<endl;
        return -1;
    }

    double minLength = 100.0;
    std::pair<Point, Point> smallestEdge = {polygon[0], polygon[0]};

    size_t n = polygon.size();
    for (size_t i = 0; i < n; ++i) {
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % n]; // Wrap around to the first point

        // Calculate the distance between p1 and p2
        double length = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));

        // Update the largest edge if this one is longer
        if (length < minLength) {
            minLength = length;
            smallestEdge = {p1, p2};
        }
    }

    return minLength;
}
*/
bool CityGen::findSmallestEdgeAmount(vector<Point> polygon, float minEdgeLength) {
    if (polygon.size() < 2) {
        // cout << "Polygon size: " << polygon.size() << endl;
        // cerr << "Polygon must have at least two points to form an edge." << endl;
        return false; // Return false if the polygon is invalid
    }

    size_t n = polygon.size();
    int edgesBelowMin = 0; // Counter for edges below the minimum length

    for (size_t i = 0; i < n; ++i) {
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % n]; // Wrap around to the first point

        // Calculate the distance between p1 and p2
        double length = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));

        if (length < minEdgeLength) {
            edgesBelowMin++;
            if (edgesBelowMin >= 2) {
                return true; // Return true if two or more edges are below the minimum
            }
        }
    }

    return false; // Return false if fewer than two edges are below the minimum
}

#pragma endregion

#pragma region Debugging
void CityGen::addDebug1(Point p)
{
    debugs1.push_back({{p.x-5.0f,p.y -5.0f}, {p.x+5.0f,p.y +5.0f}});
    debugs1.push_back({{p.x-5.0f,p.y +5.0f}, {p.x +5.0f,p.y -5.0f}});
}
void CityGen::addLineDebug1(Point p1,Point p2)
{
    debugs1.push_back({p1,p2});
}
void CityGen::addLineToVector1(Point p1,Point vector)
{
    addDebug1(p1);
    addLineDebug1(p1,{p1.x + vector.x * 10.0f, p1.y + vector.y * 10.0f});
}
void CityGen::addDebug(Point p)
{
    debugs2.push_back({{p.x-5.0f,p.y -5.0f}, {p.x+5.0f,p.y +5.0f}});
    debugs2.push_back({{p.x-5.0f,p.y +5.0f}, {p.x +5.0f,p.y -5.0f}});
}
void CityGen::addLineDebug(Point p1,Point p2)
{
    debugs2.push_back({p1,p2});
}
void CityGen::addLineToVector(Point p1,Point vector)
{
    addDebug(p1);
    addLineDebug(p1,{p1.x + vector.x * 10.0f, p1.y + vector.y * 10.0f});
}
vector<Point> CityGen::findIntersectionsWithBoundary(const Line& line) {
    vector<Point> intersections;

    // Define boundary lines
    Line top = {0, 1, -maxY};          // Horizontal: y = maxY
    Line bottom = {0, 1, -minY};        // Horizontal: y = minY
    Line left = {1, 0, -minX};         // Vertical: x = minX
    Line right = {1, 0, -maxX};         // Vertical: x   = maxX

    // Check intersection boundary
    Point intersection;
    if (findLineIntersection(line, top, intersection)) {
        if (intersection.x >= minX && intersection.x <= maxX) {
            intersections.push_back(intersection);
        }
    }
    if (findLineIntersection(line, bottom, intersection)) {
        if (intersection.x >= minX && intersection.x <= maxX) {
            intersections.push_back(intersection);
        }
    }
    if (findLineIntersection(line, left, intersection)) {
        if (intersection.y >= minY && intersection.y <= maxY) {
            intersections.push_back(intersection);
        }
    }
    if (findLineIntersection(line, right, intersection)) {
        if (intersection.y >= minY && intersection.y <= maxY) {
            intersections.push_back(intersection);
        }
    }

    return intersections;
}

void CityGen::drawCrosswalk(Point base, Point direction, float lineLength, float spacing, float lineHeight) {
    GLubyte r = 255, g = 255, b = 255, a = 255;
    Point perpendicular = perpendicularVector(direction.x, direction.y);

    for (float offset = -(spacing * 5); offset <= spacing * 5; offset += spacing) {
        Point lineStart = movePointInDirection(base, perpendicular, offset);
        Point lineEnd = movePointInDirection(lineStart, direction, lineLength);

        Vertex v1 = { static_cast<GLfloat>(lineStart.x), lineHeight, static_cast<GLfloat>(lineStart.y), r, g, b, a,
        0.0f,0.0f };
        Vertex v2 = { static_cast<GLfloat>(lineEnd.x), lineHeight, static_cast<GLfloat>(lineEnd.y), r, g, b, a,
        1.0f,1.0f };

        m_lines.push_back(v1);
        m_lines.push_back(v2);
    }
}

void CityGen::addRoadDecals(Point p1, Point p2) {
    float lineHeight = 0.01f;
    float lineLength = 5.0f;
    float gapLength = 8.0f;
    float intersection = 15.0f;
    float crosswalkWidth = 8.0f; // Width of crosswalk
    float crosswalkSpacing = 1.5f; // Spacing between lines in the crosswalk
    float crosswalkLineLength = 5.0f; // Length of each line in the crosswalk
    float minRoadLength = 30.0f;

    float totalLength = distanceBetweenPoints(p1, p2);
    if ((totalLength * 0.65) <= (intersection * 2)
    ||   totalLength <= minRoadLength) {
        return;
    }

    Point p1Dir = getDirectionVector(p1, p2);
    Point p2Dir = getDirectionVector(p2, p1);

    Point startCrosswalk = movePointInDirection(p1, p1Dir, -crosswalkWidth);
    Point endCrosswalk = movePointInDirection(p2, p2Dir, -crosswalkWidth);

    p1 = movePointInDirection(p1, p1Dir, intersection);
    p2 = movePointInDirection(p2, p2Dir, intersection);

    // Add crosswalk at the start
    drawCrosswalk(p1, p1Dir, crosswalkLineLength, crosswalkSpacing, lineHeight);

    // Add crosswalk at the end
    drawCrosswalk(p2, p2Dir, crosswalkLineLength, crosswalkSpacing, lineHeight);

    totalLength = distanceBetweenPoints(p1, p2);

    float segmentLength = lineLength + gapLength;
    int numSegments = static_cast<int>(totalLength / segmentLength);

    GLubyte r = 255, g = 255, b = 255, a = 255;

    Point direction = {
        (p2.x - p1.x) / totalLength,
        (p2.y - p1.y) / totalLength
    };

    for (int i = 0; i <= numSegments; ++i) {
        float startDist = i * segmentLength;
        float endDist = startDist + lineLength;

        if (startDist >= totalLength) break;

        // Clamp end distance to total length
        if (endDist > totalLength) endDist = totalLength;

        // Calculate start and end points of the solid line segment
        Point interpStart = {
            p1.x + startDist * direction.x,
            p1.y + startDist * direction.y
        };

        Point interpEnd = {
            p1.x + endDist * direction.x,
            p1.y + endDist * direction.y
        };

        // Add the solid segment as a pair of vertices
        Vertex v1 = { static_cast<GLfloat>(interpStart.x), lineHeight, static_cast<GLfloat>(interpStart.y), r, g, b, a,
        0.0f,0.0f };
        Vertex v2 = { static_cast<GLfloat>(interpEnd.x), lineHeight, static_cast<GLfloat>(interpEnd.y), r, g, b, a,
        1.0f,1.0f };

        m_lines.push_back(v1);
        m_lines.push_back(v2);
    }
}

/*void CityGen::addRoadDecals(Point p1, Point p2) {
    float lineHeight = 0.01f;
    float lineLength = 5.0f;
    float gapLength = 8.0f;
    float intersection = 20.0f;

    
    float totalLength = distanceBetweenPoints(p1, p2);
    if(totalLength <= intersection * 2) {return;}

    Point p1Dir = getDirectionVector(p1,p2);
    Point p2Dir = getDirectionVector(p2,p1);

    p1 = movePointInDirection(p1,p1Dir,intersection);
    p2 = movePointInDirection(p2,p2Dir,intersection);

    totalLength = distanceBetweenPoints(p1, p2);

    float segmentLength = lineLength + gapLength;
    int numSegments = static_cast<int>(totalLength / segmentLength);

    GLubyte r = 255, g = 255, b = 255, a = 255;

    Point direction = {
        (p2.x - p1.x) / totalLength,
        (p2.y - p1.y) / totalLength
    };

    for (int i = 0; i <= numSegments; ++i) {
        float startDist = i * segmentLength;
        float endDist = startDist + lineLength;

        if (startDist >= totalLength) break;

        // Clamp end distance to total length
        if (endDist > totalLength) endDist = totalLength;

        // Calculate start and end points of the solid line segment
        Point interpStart = {
            p1.x + startDist * direction.x,
            p1.y + startDist * direction.y
        };

        Point interpEnd = {
            p1.x + endDist * direction.x,
            p1.y + endDist * direction.y
        };

        // Add the solid segment as a pair of vertices
        Vertex v1 = { static_cast<GLfloat>(interpStart.x), lineHeight, static_cast<GLfloat>(interpStart.y), r, g, b, a };
        Vertex v2 = { static_cast<GLfloat>(interpEnd.x), lineHeight, static_cast<GLfloat>(interpEnd.y), r, g, b, a };

        m_lines.push_back(v1);
        m_lines.push_back(v2);
    }
}
*/
#pragma endregion

#pragma region Voronoi gen
vector<Point> CityGen::generateSites(int numSites, unsigned int seed) {

    mt19937 gen(seed);
    uniform_real_distribution<double> distX(minX + 2, maxX - 2);
    uniform_real_distribution<double> distY(minY + 2, maxY - 2);

    vector<Point> sites;
    for (int i = 0; i < numSites; ++i) {
        float x = distX(gen);
        float y = distY(gen);
        sites.push_back({x, y});
    }
    return sites;
}
void CityGen::computeVoronoiDiagram() {
    m_voronoiCells.clear();
    m_voronoiCells.resize(m_sites.size());

    vector<Point> boundingPolygon = {
        { minX, minY },
        { maxX, minY },
        { maxX, maxY },
        { minX, maxY }
    };

    for (size_t i = 0; i < m_sites.size(); i++) {
        const Point& site = m_sites[i];
        vector<Point> cell = boundingPolygon;

        for (size_t j = 0; j < m_sites.size(); j++) {
            if (i == j) continue;

            const Point& otherSite = m_sites[j];

            Line bisector = perpendicularBisector(site, otherSite);

            double eval = evaluate(bisector, site);
            bool keepPositiveSide = eval > 0;
            
            // cell = clipPolygon(cell, bisector);
            keepPositiveSide ? cell = clipPolygon(cell,bisector).first : cell = clipPolygon(cell,bisector).second; //[0] positive, [1] negative

            if (cell.empty()) {
                cout<<"\nCell is empty somefuckinghow";
                break;
            }
        }
        m_voronoiCells[i] = cell;
    }
}
void CityGen::computeChunks() {
    m_chunks = m_voronoiCells;

    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    mt19937 gen(seed);
    std::uniform_int_distribution<int> chanceDist(1, 8);

    float minimumArea = 600.0f;
    float minimumEdgeLength = 50.0f;

    size_t currentChunkIndex = 0;
    vector<vector<Point>> tempPolygonList;

    for(auto currentChunk : m_chunks) {
        // vector<Point> currentChunk = m_chunks[currentChunkIndex];
        if (calculatePolygonArea(currentChunk) <= minimumArea
         || currentChunk.size() <= 3)
        {
            tempPolygonList.push_back(currentChunk);
            continue;
        }

        pair<Point, Point> largestEdge = findLargestEdge(currentChunk);
        float edgeLength = distanceBetweenPoints(largestEdge.first, largestEdge.second);

        // Skip if the largest edge is too short
        if (edgeLength <= minimumEdgeLength) {
            tempPolygonList.push_back(currentChunk);
            continue;
        }

        uniform_int_distribution<int> splitDist(0, static_cast<int>(currentChunk.size() - 1));
        int index = splitDist(gen);
        // int index = 1; //removed randomizing
        int nextIndex = (index + 1) % currentChunk.size();

        Line cut = CreateLineFromPoints(currentChunk[index], currentChunk[nextIndex]);
        cut = moveLineToCenter(cut, currentChunk);

        auto clipped = clipPolygon(currentChunk, cut);
        // vector<Point> positive = clipped.first;
        // vector<Point> negative = clipped.second;
        tempPolygonList.push_back(clipped.first);
        tempPolygonList.push_back(clipped.second);

        // if (!positive.empty()) {
        //     tempPolygonList.push_back(positive);
        // }
        // if (!negative.empty()) {
        //     tempPolygonList.push_back(negative);
        // }
    }
    m_chunks = tempPolygonList;

    m_voronoiCells = tempPolygonList;

    // for (size_t i = 0; i < m_chunks.size(); i++) {
    //     vector<Point> x = m_chunks[i];
    //     m_chunks[i] = scalePolygon(x, 0.85f);
    //     // m_blocks[i] = scalePolygon(x, 0.85f);
    // }
    /*
    m_chunks = m_voronoiCells;

    float minimumArea = 600.0f;
    float minimumEdgeLength = 20.0f;

    size_t currentChunkIndex = 0;

    while (currentChunkIndex < m_chunks.size()) {
        vector<Point> currentChunk = m_chunks[currentChunkIndex];

        // Skip small polygons
        if (calculatePolygonArea(currentChunk) <= minimumArea || currentChunk.size() < 3) {
            currentChunkIndex++;
            continue;
        }

        // Find the largest edge of the polygon
        pair<Point, Point> largestEdge = findLargestEdge(currentChunk);
        float edgeLength = distanceBetweenPoints(largestEdge.first, largestEdge.second);

        // Skip if the largest edge is too short
        if (edgeLength <= minimumEdgeLength) {
            currentChunkIndex++;
            continue;
        }

        // Create a perpendicular bisector for the largest edge
        Line bisector = perpendicularBisector(largestEdge.first, largestEdge.second);
        bisector = moveLineToCenter(bisector, currentChunk);

        // Split the polygon using the bisector
        auto clipped = clipPolygon(currentChunk, bisector);
        vector<Point> positive = clipped.first;
        vector<Point> negative = clipped.second;

        // Add the resulting polygons back if they're valid
        if (!positive.empty()) {
            m_chunks.push_back(positive);
        }
        if (!negative.empty()) {
            m_chunks.push_back(negative);
        }

        // Remove the current chunk
        m_chunks.erase(m_chunks.begin() + currentChunkIndex);
    }*/
    /*
    // m_voronoiCells = m_chunks;
    m_chunks = m_voronoiCells;

    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    mt19937 gen(seed);
    std::uniform_int_distribution<int> chanceDist(1, 8);

    float minimumArea = 600.0f;

    size_t currentChunkIndex = 0;

    while (currentChunkIndex < m_chunks.size()) {
        vector<Point> currentChunk = m_chunks[currentChunkIndex];
        if (calculatePolygonArea(currentChunk) <= minimumArea)
        {
            currentChunkIndex++; continue;
        }
        
        if (currentChunk.empty() || (chanceDist(gen) == 1)) // 1 in 8 chance of skipping step
        {currentChunkIndex++; continue;}

        if (calculatePolygonArea(currentChunk) >= maximumChunkSize) {
            // cout << "Splitting chunk at index: " << currentChunkIndex << endl;
            if (currentChunk.size() <= 3) {
                // cout<<"Chunk too small to cut"<<endl;
                currentChunkIndex++;
                continue;
            }

            uniform_int_distribution<int> splitDist(0, static_cast<int>(currentChunk.size() - 1));

            // int index = splitDist(gen);
            int index = 1; //removed randomizing
            int nextIndex = (index + 1) % currentChunk.size();

            Line cut = CreateLineFromPoints(currentChunk[index], currentChunk[nextIndex]);
            cut = moveLineToCenter(cut, currentChunk);

            auto clipped = clipPolygon(currentChunk, cut);
            vector<Point> positive = clipped.first;
            vector<Point> negative = clipped.second;

            if (!positive.empty()) {
                m_chunks.push_back(positive);
            }
            if (!negative.empty()) {
                m_chunks.push_back(negative);
            }

            m_chunks.erase(m_chunks.begin() + currentChunkIndex);
        }
        else {
            // Current chunk is within the size limit; move to the next
            currentChunkIndex++;
        }
    }

    // m_voronoiCells = m_chunks;

    // for (size_t i = 0; i < m_chunks.size(); i++) {
    //     vector<Point> x = m_chunks[i];
    //     m_chunks[i] = scalePolygon(x, 0.85f);
    //     // m_blocks[i] = scalePolygon(x, 0.85f);
    // }

    */
}
void CityGen::sweepToBlocks()
{
    m_blocks.clear();
    m_blocks.resize(m_chunks.size());
    m_buildings.clear();
    
    m_chunks = m_voronoiCells;

    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];
        m_chunks[i] = scalePolygon(x, 0.85f);
        // m_blocks[i] = scalePolygon(x, 0.85f);
    }
    
  // Cuts chunks into blocks
    for (size_t i = 0; i < m_chunks.size(); i++) { //m_chunks.size()
        #pragma region prep
        vector<Point> x = m_chunks[i];

        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(minMoveAmount, maxMoveAmount); // Range for moveAmount

        pair<Point,Point> largestPointPair = findLargestEdge(x);
        
        double edgeLength = distanceBetweenPoints(largestPointPair.first,largestPointPair.second);

        Line largestEdge = CreateLineFromPoints(largestPointPair.first,largestPointPair.second);

        Point midpoint = findMidpoint(largestPointPair.first,largestPointPair.second);
        
        Point centerPoint = getCentroid(x);
        
        // Point directionVector = getDirectionVector(midpoint,centerPoint);
        Point directionVector = getPerpendicularDirVector(largestPointPair,centerPoint);

        Point goToPoint = midpoint;
        largestEdge = moveToPoint(largestEdge,goToPoint);

        int iterationAmount = 0;

        double eval = evaluate(largestEdge, centerPoint);
        bool keepPositiveSide = eval > 0;
        
        m_strips.clear();
        #pragma endregion
        #pragma region Part1
        while(true)
        {
            iterationAmount++;

            moveAmount = dis(gen);

            goToPoint = movePointInDirection(goToPoint, directionVector, moveAmount);
            largestEdge = moveToPoint(largestEdge,goToPoint);
            auto clipped = clipPolygon(x, largestEdge);
            vector<Point> positive = clipped.first;
            vector<Point> negative = clipped.second;

            if((positive.empty() && negative.empty())) {break;}

            if((findSmallestEdgeAmount(positive, minEdge) && positive.size() <= 4) 
            ||(findSmallestEdgeAmount(negative, minEdge) && negative.size() <= 4))
            {
                m_strips.push_back(x);
                break;
            }
            else if(keepPositiveSide)
            {
                x = positive;

                if(!negative.empty()) {m_strips.push_back(negative);}//m_blocks.push_back(negative);}
            }else
            {
                x = negative;
                //m_strips
                if(!positive.empty()) {m_strips.push_back(positive);}//{m_blocks.push_back(positive);}
            }
        }
        
        //You have the building strips... Now cut through them.
        //m_strips
        #pragma endregion
        #pragma region prep
        
        float distance = edgeLength * .7f;

        iterationAmount = 0;

        // Create a perpendicular direction vector
        directionVector = perpendicularVector(directionVector.x, directionVector.y);
        Point oppVector = oppositeVector(directionVector.x, directionVector.y);

        goToPoint = centerPoint;
        // addLineToVector(goToPoint,oppVector);
        goToPoint = movePointInDirection(goToPoint, oppVector, distance);
        // addLineToVector(goToPoint,directionVector);
        
        std::cout<<"\n\nDistance: "<<distance<<endl;

        // largestEdge = moveToPoint(largestEdge, centerPoint);
        largestEdge = makePerpendicularLine(largestEdge);
        largestEdge = moveToPoint(largestEdge, goToPoint);
        // debugs1.push_back(findIntersectionsWithBoundary(largestEdge));

        eval = evaluate(largestEdge, centerPoint);
        keepPositiveSide = eval > 0;

        std::cout << "\n\nDebug: START" << endl;

        #pragma endregion
        #pragma region Part2

        vector<vector<Point>> tempStrips;
        
        while(!m_strips.empty())
        {
            iterationAmount++;
            // cout<<"\n\n\nSTRIP SIZE: "<<m_strips.size()<<"\n\n"<<endl;

            moveAmount = dis(gen);

            goToPoint = movePointInDirection(goToPoint, directionVector, moveAmount);
            // addLineToVector(goToPoint,directionVector);
            largestEdge = moveToPoint(largestEdge, goToPoint);

            // debugs1.push_back(findIntersectionsWithBoundary(largestEdge));

            for (size_t j = 0; j < m_strips.size(); j++)
            {
                // cout<<"Loop: "<<j<<"\n"<<endl;

                vector<Point> currPolygon = m_strips[j];
    
                auto clipped = clipPolygon(currPolygon, largestEdge);
                vector<Point> positive = clipped.first;
                vector<Point> negative = clipped.second;
    
                // cout << " - Positive side: " << positive.size() <<
                // " , Negative side: " << negative.size() << endl;
    
                if((positive.empty() && negative.empty())) {
                    // cout << "Both empty" << j + 1 << ". Skipping." << endl;
                    continue;
                }
    
                bool positiveMeetsCriteria = (findSmallestEdgeAmount(positive, minEdge2) && positive.size() <= 4);
                bool negativeMeetsCriteria = (findSmallestEdgeAmount(negative, minEdge2) && negative.size() <= 4);
                // bool positiveMeetsCriteria = (calculatePolygonArea(positive) <= minPolygonArea ||
                // (findSmallestEdgeAmount(negative, minEdge2) && positive.size() <= 4));

                // bool negativeMeetsCriteria = (calculatePolygonArea(negative) <= minPolygonArea ||
                // (findSmallestEdgeAmount(negative, minEdge2) && negative.size() <= 4));

                // cout << "Positive: " << (positiveMeetsCriteria ? "True" : "False") 
                //      << ", Negative: " << (negativeMeetsCriteria ? "True" : "False") << endl;
                    
                if((positiveMeetsCriteria) || (negativeMeetsCriteria))
                {
                    if(keepPositiveSide && !positive.empty())
                    {
                        currPolygon = positive;
                        tempStrips.push_back(currPolygon);
                    }
                    else if(!keepPositiveSide && !negative.empty())
                    {
                        currPolygon = negative;
                        tempStrips.push_back(currPolygon);
                    }
                    else
                    {
                        m_buildings.push_back(currPolygon);
                        // cout << "Criteria MET Added strip " << j + 1 << " as a building." << endl;
                    }
                }
                else if(keepPositiveSide)
                {
                    if(positive.empty()) {
                        // cout << "Positive empty for strip " << j + 1 << ". Skipping." << endl;
                    }else
                    {
                        currPolygon = positive;
                        tempStrips.push_back(currPolygon);
                    }

                    if(!negative.empty()) { //if building clipped, add to buildings
                        if(calculatePolygonArea(negative) > minPolygonArea)
                        {
                            m_buildings.push_back(negative);
                            // cout << "negative side as a building for strip " << j + 1 << "." << endl;

                        }                        
                    }
                }
                else
                {
                    if(negative.empty()) {
                        // cout << "Negative empty for strip " << j + 1 << ". Skipping." << endl;
                    }
                    else
                    {
                        currPolygon = negative;
                        tempStrips.push_back(currPolygon);
                    }

                    if(!positive.empty()){
                        if(calculatePolygonArea(positive) > minPolygonArea)
                        {
                            m_buildings.push_back(positive);
                            // cout << "positive side as a building for strip " << j + 1 << "." << endl;
                        }
                    }

                }
            }
            m_strips = tempStrips;
            tempStrips.clear();
        }
        #pragma endregion
    }
}
#pragma endregion

#pragma endregion

#pragma region Render

void CityGen::pushVertexData(GLuint vao, GLuint vbo, vector<Vertex> vertices)
{
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, x));                      // Position
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (void*)offsetof(Vertex, r));
    glEnableVertexAttribArray(1);
    
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, u));
    glEnableVertexAttribArray(2);

    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, nx));
    glEnableVertexAttribArray(3);
    
    glBindVertexArray(0);
}

void CityGen::buildVertexData() {
    // /*
    m_vertices.clear();
    // m_lines.clear();
/*
    for (size_t i = 0; i < m_buildings.size(); i++) {
        vector<Point> x = m_buildings[i];
    
        m_buildings[i] = scalePolygon(x, 0.85f);
        //write new polygon scaling logic;
        //Have each vertex move towards centroid by a certain amount,
        //if amount >= distance between centroid and vertex cancel
        //go through all verticies and delete ones that are close to eachother
        //Can store vertex before and compare after.
    }
    
    vector<GLubyte> colors = {
    255, 0, 0,     // Bright Red
    0, 255, 0,     // Bright Green
    0, 0, 255,     // Bright Blue
    255, 255, 0,   // Bright Yellow
    255, 0, 255,   // Bright Magenta
    0, 255, 255,   // Bright Cyan
    191, 0, 191,   // Vibrant Purple
    255, 140, 0,   // Vibrant Orange
    0, 200, 0,     // Bright Green (adjusted)
    200, 200, 200  // Light Gray (stands out against black)
    };

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> heightDist(30.0f, 100.0f); 
    std::uniform_int_distribution<int> chanceDist(1, 30);
    // float wallHeight = 40.0f;
    
    vector<Vertex> topVertices;
    for (size_t i = 0; i < m_buildings.size(); ++i) { //m_buildings.size()
        topVertices.clear();
        const vector<Point>& cell = m_buildings[i];
        
        float wallHeight = (chanceDist(gen) == 1) ? 120.0f : heightDist(gen);

        
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        
        GLubyte a = 255;

        float texWidth = 10.0f; 
        float texHeight = 10.0f;
        
        //WALL VERTS
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];
            float width = distanceBetweenPoints(p1,p2);

            float tileU = width / texWidth;
            float tileV = wallHeight / texHeight;

            // tileU = 1.0f;
            // tileV = 1.0f;

            Vertex v1 = { static_cast<GLfloat>(p1.x), 0.0f, static_cast<GLfloat>(p1.y), r, g, b, a,
            0.0f,0.0f };
            Vertex v2 = { static_cast<GLfloat>(p2.x), 0.0f, static_cast<GLfloat>(p2.y), r, g, b, a,
            tileU,0.0f };

            Vertex top1 = { static_cast<GLfloat>(p1.x), wallHeight, static_cast<GLfloat>(p1.y), r, g, b, a,
            0.0f,tileV };
            Vertex top2 = { static_cast<GLfloat>(p2.x), wallHeight, static_cast<GLfloat>(p2.y), r, g, b, a,
            tileU,tileV};

            // Vertex v1 = { static_cast<GLfloat>(p1.x), 0.0f, static_cast<GLfloat>(p1.y), r, g, b, a, 
            // 0.0f,0.0f };
            // Vertex v2 = { static_cast<GLfloat>(p2.x), 0.0f, static_cast<GLfloat>(p2.y), r, g, b, a,
            // 1.0f,0.0f };
            // Vertex top1 = { static_cast<GLfloat>(p1.x), wallHeight, static_cast<GLfloat>(p1.y), r, g, b, a,
            // 0.0f,1.0f };
            // Vertex top2 = { static_cast<GLfloat>(p2.x), wallHeight, static_cast<GLfloat>(p2.y), r, g, b, a,
            // 1.0f,1.0f };

            m_vertices.push_back(v1);
            m_vertices.push_back(top1);
            m_vertices.push_back(v2);
            // m_vertices.push_back(top1);

            m_vertices.push_back(top1);
            m_vertices.push_back(top2);
            m_vertices.push_back(v2);
            // m_vertices.push_back(top2);
            
            // top1.y += 0.1f;
            topVertices.push_back(top1);
        }
        if(topVertices.empty()){continue;}
        //ROOF - TOP VERTS
        for (size_t j = 1; j < topVertices.size() - 1; ++j) {
            Vertex v1 = topVertices[0];
            Vertex v2 = topVertices[j];
            Vertex v3 = topVertices[j + 1];
            // Add roof triangles
            m_vertices.push_back(v1);
            m_vertices.push_back(v3);
            m_vertices.push_back(v2);
            // m_vertices.push_back(v3);
        }
    }
    
    m_lines.clear();
    for (size_t i = 0; i < m_voronoiCells.size(); ++i) 
    {
        const vector<Point>& cell = m_voronoiCells[i];     
        GLubyte r = 255;
        GLubyte g = 255;
        GLubyte b = 255;
        GLubyte a = 255;
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];
            // Vertex v1 = { static_cast<GLfloat>(p1.x), -20.0f, static_cast<GLfloat>(p1.y), r,g,b,a,0.0f,0.0f};
            // Vertex v2 = { static_cast<GLfloat>(p2.x), -20.0f, static_cast<GLfloat>(p2.y), r,g,b,a,0.0f,0.0f};
            // m_lines.push_back(v1);
            // m_lines.push_back(v2);
            
            addRoadDecals(p1,p2);
        }
    }
*/
    if(Debug)
    {   
        // for (size_t j = 0; j < m_voronoiCells.size(); ++j) //debugs
        // {
        //     if(m_voronoiCells[j].empty()) continue;
        //     Point p1 = m_voronoiCells[j][0];
        //     Point p2 = m_voronoiCells[j][1];
        //     Vertex v1 = { static_cast<GLfloat>(p1.x), 10.0f, static_cast<GLfloat>(p1.y), 255, 0, 0, 255};
        //     Vertex v2 = { static_cast<GLfloat>(p2.x), 10.0f, static_cast<GLfloat>(p2.y), 255, 0, 0, 255};
        //     m_lines.push_back(v1);
        //     m_lines.push_back(v2);
        // }

        /*
       for (size_t i = 0; i < m_voronoiCells.size(); ++i) 
        {
            const vector<Point>& cell = m_voronoiCells[i];     
            GLubyte r = 255;
            GLubyte g = 255;
            GLubyte b = 255;
            GLubyte a = 255;
            for (size_t j = 0; j < cell.size(); ++j) {
                Point p1 = cell[j];
                Point p2 = cell[(j + 1) % cell.size()];
                addLineDebug1(p1,p2);
                // addRoadDecals(p1,p2);
            }
        }*/
        // for (size_t j = 0; j < debugs1.size(); ++j)
        // {
        //     if(debugs1[j].empty()) continue;
        //     Point p1 = debugs1[j][0];
        //     Point p2 = debugs1[j][1];
        //     Vertex v1 = { static_cast<GLfloat>(p1.x), 20.0f, static_cast<GLfloat>(p1.y), 0, 0, 255, 255};
        //     Vertex v2 = { static_cast<GLfloat>(p2.x), 20.0f, static_cast<GLfloat>(p2.y), 0, 0, 255, 255};
        //     // m_lines.push_back(v1);
        //     // m_lines.push_back(v2);
        // }
        /*
        for (size_t j = 0; j < debugs2.size(); ++j)
        {
            if(debugs2[j].empty()) continue;
            Point p1 = debugs2[j][0];
            Point p2 = debugs2[j][1];
            Vertex v1 = { static_cast<GLfloat>(p1.x), 30.0f, static_cast<GLfloat>(p1.y), 0, 255, 0, 255};
            Vertex v2 = { static_cast<GLfloat>(p2.x), 30.0f, static_cast<GLfloat>(p2.y), 0, 255, 0, 255};
            m_lines.push_back(v1);
            m_lines.push_back(v2);
        }
        */
    }

// /*UV DEBUG

    vec3 normal = vec3(-1.0f,0.0f,-1.0f);

    Vertex x00 = { static_cast<GLfloat>(-20.0f), 20.0f, static_cast<GLfloat>(0.0f),
    255, 0, 0, 255,
    0.0f,0.0f
    ,normal.x,normal.y,normal.z};

    Vertex x10 = { static_cast<GLfloat>(20.0f), 20.0f, static_cast<GLfloat>(0.0f),
    255, 0, 0, 255,
    1.0f,0.0f
    ,normal.x,normal.y,normal.z};

    Vertex x01 = { static_cast<GLfloat>(-20.0f), 60.0f, static_cast<GLfloat>(0.0f),
    255, 0, 0, 255,
    0.0f,1.0f
    ,normal.x,normal.y,normal.z};

    Vertex x11 = { static_cast<GLfloat>(20.0f), 60.0f, static_cast<GLfloat>(0.0f),
    255, 0, 0, 255,
    1.0f,1.0f
    ,normal.x,normal.y,normal.z};

    m_vertices.push_back(x00);
    m_vertices.push_back(x10);
    m_vertices.push_back(x11);
    
    m_vertices.push_back(x00);
    m_vertices.push_back(x11);
    m_vertices.push_back(x01);

    std::cout << "Size of Vertex: " << sizeof(Vertex) << " bytes" << std::endl;
    std::cout << "Offset of r: " << offsetof(Vertex, r) << " bytes" << std::endl;
    std::cout << "Offset of u: " << offsetof(Vertex, u) << " bytes" << std::endl;
    std::cout << "Offset of nx: " << offsetof(Vertex, nx) << " bytes" << std::endl;

// */

    #pragma region Vertex Prep
    m_numVertices = static_cast<int>(m_vertices.size());// * 5.0f); // sussyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

    // Create and bind VAO and VBO
    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1, &m_vbo);
    pushVertexData(m_vao,m_vbo,m_vertices);
    
    //Debugging UV
    // GLint positionLoc = glGetAttribLocation(m_program, "a_position");
    // GLint texcoordLoc = glGetAttribLocation(m_program, "a_texcoord");
    // std::cout << "Position Attribute Location: " << positionLoc << std::endl;
    // std::cout << "Texcoord Attribute Location: " << texcoordLoc << std::endl;
    std::cout << "sizeof(Vertex): " << sizeof(Vertex) << std::endl;
    std::cout << "Offset of Position: " << offsetof(Vertex, x) << std::endl;
    std::cout << "Offset of UV: " << offsetof(Vertex, u) << std::endl;
    std::cout << "Buffer Data Size: " << m_vertices.size() * sizeof(Vertex) << " bytes" << std::endl;
    std::cout << "Stride: " << sizeof(Vertex) << " bytes" << std::endl;
    //---

    // LINES - Unbind VAO
    glBindVertexArray(0);
    m_numLines = static_cast<int>(m_lines.size());
    glGenVertexArrays(1, &m_vaoLines);
    glGenBuffers(1, &m_vboLines);
    pushVertexData(m_vaoLines,m_vboLines,m_lines);

    glBindVertexArray(0);

    m_buildingTexture = wolf::TextureManager::CreateTexture("data/building1.tga");
    m_buildingTexture->SetWrapMode(wolf::Texture::WrapMode::WM_Repeat,wolf::Texture::WrapMode::WM_Repeat);
    m_buildingTexture->SetFilterMode(wolf::Texture::FilterMode::FM_LinearMipmap, wolf::Texture::FilterMode::FM_LinearMipmap);
    // m_buildingTexture->SetFilterMode(wolf::Texture::FilterMode::FM_Linear, wolf::Texture::FilterMode::FM_Linear);

    // */
    #pragma endregion
}
void CityGen::generate(GLuint program) {
    m_program = wolf::ProgramManager::CreateProgram("data/cube.vsh", "data/cube.fsh");
    // m_program = program;
    
    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    m_sites = generateSites(numSites, seed);
    
    // m_sites = {
    // {-50,-50},
    // {50,50}};

    computeVoronoiDiagram();
    // computeChunks();
    sweepToBlocks();
    buildVertexData();
}

void CityGen::deGenerate() {
    // Clean up OpenGL resources
    if (m_program) {
        wolf::ProgramManager::DestroyProgram(m_program);
        m_program = nullptr;
    }
    if (m_buildingTexture) {
        wolf::TextureManager::DestroyTexture(m_buildingTexture);
        m_buildingTexture = nullptr;
    }
    if (m_vbo) {
        glDeleteBuffers(1, &m_vbo);
        m_vbo = 0;
    }
    if (m_vao) {
        glDeleteVertexArrays(1, &m_vao);
        m_vao = 0;
    }
    m_vertices.clear();
    debugs.clear();
    debugs1.clear();
    debugs2.clear();
    m_sites.clear();
    m_voronoiCells.clear();
    m_buildings.clear();
}
void CityGen::reGenerate() {
    deGenerate();
    generate(0);
}
void CityGen::render(const glm::mat4& proj, const glm::mat4& view, float m_timer) {
    if (!m_program || !m_buildingTexture) {
        cerr << "Shader program or texture not initialized." << endl;
        return;
    }
    m_program->Bind();

    // Set uniform variables using Wolf's SetUniform
    m_program->SetUniform("projection", proj);
    m_program->SetUniform("view", view);
    m_program->SetUniform("world", glm::mat4(1.0f)); // Assuming no world transformation
    m_program->SetUniform("u_time", m_timer);
    m_program->SetUniform("u_texture", 0); // Texture unit 0

    // Bind texture to texture unit 0
    m_buildingTexture->Bind(0);

    // Render buildings
    glBindVertexArray(m_vao);
    glDrawArrays(GL_TRIANGLES, 0, m_numVertices);
    glBindVertexArray(0);

    // Render roads (lines)
    glLineWidth(5.0f);
    glEnable(GL_LINE_SMOOTH);
    glBindVertexArray(m_vaoLines);
    glDrawArrays(GL_LINES, 0, m_numLines);
    glBindVertexArray(0);

    // Unbind the shader program
    // m_program->Bind(0);
    glUseProgram(0);
}

#pragma endregion
