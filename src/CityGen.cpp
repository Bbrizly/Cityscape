#include "CityGen.h"
#include "Building.h"
#include "Road.h"
#include "Cube.h"
#include "Vertex.h"
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <iostream>
#include <glm/glm.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
using namespace std;
using namespace glm;
vector<Cube> cubes;

#pragma region mathematical utility
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
//---

//---

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

Line moveLinePerpendicularly(const Line& line, double distance, const Point& fixedDirection) {
    // Compute the magnitude of the line's normal vector (a, b)
    double magnitude = sqrt(line.a * line.a + line.b * line.b);
    if (magnitude < 1e-9) {
        cerr << "Invalid line: Normal vector has zero length." << endl;
        return line;
    }

    // Normalize the fixed direction vector
    double dirMagnitude = sqrt(fixedDirection.x * fixedDirection.x + fixedDirection.y * fixedDirection.y);
    if (dirMagnitude < 1e-9) {
        cerr << "Invalid fixed direction: Zero-length vector." << endl;
        return line;
    }
    Point normalizedDirection = {fixedDirection.x / dirMagnitude, fixedDirection.y / dirMagnitude};

    // Project the normalized direction onto the line's normal vector
    double dotProduct = normalizedDirection.x * line.a + normalizedDirection.y * line.b;

    // Calculate signed distance
    double signedDistance = distance * (dotProduct > 0 ? 1 : -1);

    // Adjust the line's c value
    double newC = line.c + signedDistance * magnitude;

    return {line.a, line.b, newC};
}
// Line moveLinePerpendicularly(const Line& line, double distance) {
//     double magnitude = sqrt(line.a * line.a + line.b * line.b);
//     double newC = line.c + distance * magnitude; // Shift the line by the given distance
//     return {line.a, line.b, newC};
// }
#pragma endregion

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
        double length = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));

        // Update the largest edge if this one is longer
        if (length > maxLength) {
            maxLength = length;
            largestEdge = {p1, p2};
        }
    }
    
    debugs2.push_back({largestEdge.first,largestEdge.second});
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
    cout << "=1=-FUNKY Line: " << l.a << "x + " << l.b << "y + " << l.c << " = 0" << endl;
    l.c = -(l.a * p.x + l.b * p.y);
    cout << "=2=-FUNKY Line: " << l.a << "x + " << l.b << "y + " << l.c << " = 0" << endl;
    return l;
}

#pragma region clutter

#pragma region Debug
void CityGen::CreateCubesAlongLine(const Point& start, const Point& end, int height) {
    glm::vec3 direction(end.x - start.x, 0.0f, end.y - start.y);

    glm::vec3 normalizedDir = normalize(direction);

    float lineLength = length(direction);
    float spacing = 2.5f;

    for (int i = 0; i < lineLength/spacing; ++i) {
        vec3 position = vec3(start.x, 0.0f, start.y) + normalizedDir * (spacing * i);
        
        Cube c(position, vec3(1.0f,1.0f,height), vec3(0),vec4(1));
        c.init(m_program);
        cubes.push_back(c);
    }
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
#pragma endregion

#pragma region Voronoi steps
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

pair<vector<Point>, vector<Point>> CityGen::clipPolygon(vector<Point> polygon, Line l) 
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

/*
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
*/
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
#pragma endregion

vector<Line> CityGen::generateSweepLines(const Point& origin, const Point& direction, double minSpacing, double maxSpacing, double maxDistance) {
    std::vector<Line> sweepLines;
    double distance = 0.0;
    std::mt19937 gen(static_cast<unsigned int>(time(nullptr)));
    std::uniform_real_distribution<double> distSpacing(minSpacing, maxSpacing);
    
    // Generate sweep lines until maxDistance is reached
    while (distance < maxDistance) {
        distance += distSpacing(gen);
        
        // Calculate a point along the perpendicular to the sweep direction
        Point perp = {-(direction.y), direction.x};
        Point sweepPoint = { origin.x + perp.x * distance, origin.y + perp.y * distance };
        
        // Define the sweep line in standard form (a*x + b*y + c = 0)
        // The sweep line is perpendicular to the sweep direction
        Line sweepLine;
        sweepLine.a = perp.x;
        sweepLine.b = perp.y;
        sweepLine.c = -(sweepLine.a * sweepPoint.x + sweepLine.b * sweepPoint.y);
        
        sweepLines.push_back(sweepLine);
    }
    
    return sweepLines;
}
vector<vector<Point>> CityGen::splitToBlocks(vector<Point>& polygon, vector<Line>& sweepLines)
{
    vector<vector<Point>> blocks;
    vector<vector<Point>> currentPolygons = { polygon };
    
    for (const auto& sweepLine : sweepLines) {
        vector<vector<Point>> newBlocks;
        for (const auto& currentPolygon : currentPolygons) {
            
            auto clipped = clipPolygon(currentPolygon, sweepLine);
            vector<Point> positive = clipped.first;
            vector<Point> negative = clipped.second;
            
            if (!positive.empty()) {
                newBlocks.push_back(positive);
            }
            if (!negative.empty()) {
                newBlocks.push_back(negative);
            }
        }
        currentPolygons = newBlocks;
    }
    
    blocks = currentPolygons;
    return blocks;
}

#pragma region checks
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
Point CityGen::normalizeVector(double x, double y) {
    double mag = std::sqrt(x * x + y * y);
    if (mag == 0.0) return {0.0, 0.0};
    return {x / mag, y / mag};
}
Point CityGen::perpendicularVector(double x, double y) {
    return {-y, x};
}
#pragma endregion

vector<Point> CityGen::offsetPolygonInward(const vector<Point>& polygon, double offsetDistance) {
    if (polygon.size() < 3) return {};

    std::vector<Point> offsetPolygon;
    size_t n = polygon.size();

    for (size_t i = 0; i < n; ++i) {
        const Point& prev = polygon[(i + n - 1) % n];
        const Point& curr = polygon[i];
        const Point& next = polygon[(i + 1) % n];
        
        // Calculate the normals for the previous and next edges
        double dx1 = curr.x - prev.x;
        double dy1 = curr.y - prev.y;
        Point normal1 = perpendicularVector(dx1, dy1);
        normal1 = normalizeVector(normal1.x, normal1.y);
        
        double dx2 = next.x - curr.x;
        double dy2 = next.y - curr.y;
        Point normal2 = perpendicularVector(dx2, dy2);
        normal2 = normalizeVector(normal2.x, normal2.y);
        
        // Average the normals
        Point avgNormal = { (normal1.x + normal2.x) / 2.0, (normal1.y + normal2.y) / 2.0 };
        avgNormal = normalizeVector(avgNormal.x, avgNormal.y);
        
        // Move the current point inward along the average normal
        Point offsetPoint = { curr.x + avgNormal.x * offsetDistance, curr.y + avgNormal.y * offsetDistance };
        offsetPolygon.push_back(offsetPoint);
    }

    // Check if the offset polygon is still valid and convex
    if (!isConvex(offsetPolygon)) {
        cout<<"HEY CONVEX ISSUE"<<endl;
        offsetPolygonInward(polygon, offsetDistance/2.0f);
        // return polygon;
    }

    return offsetPolygon;
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

Line CityGen::moveLineToCenter(Line l, vector<Point> polygon)
{
    Point centroid = getCentroid(polygon);
    return (moveToPoint(l, centroid));
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
        cout << "Polygon size: " << polygon.size() << endl;
        cerr << "Polygon must have at least two points to form an edge." << endl;
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

void CityGen::addDebug1(Point p)
{
    debugs1.push_back({{p.x-5.0f,p.y -5.0f}, {p.x+5.0f,p.y +5.0f}});
    debugs1.push_back({{p.x-5.0f,p.y +5.0f}, {p.x +5.0f,p.y -5.0f}});
}

void CityGen::addLineDebug1(Point p1,Point p2)
{
    debugs1.push_back({p1,p2});
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

void CityGen::computeChunks() {
    m_chunks = m_voronoiCells;

    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    mt19937 gen(seed);
    std::uniform_int_distribution<int> chanceDist(1, 8);

    size_t currentChunkIndex = 0;

    while (currentChunkIndex < m_chunks.size()) {
        vector<Point> currentChunk = m_chunks[currentChunkIndex];
        
        if (currentChunk.empty() || (chanceDist(gen) == 1)) // 1 in 8 chance of skipping step
        {currentChunkIndex++; continue;}

        if (calculatePolygonArea(currentChunk) >= maximumChunkSize) {
            cout << "Splitting chunk at index: " << currentChunkIndex << endl;
            if (currentChunk.size() <= 3) {
                cout<<"Chunk too small to cut"<<endl;
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

    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];
    
        m_chunks[i] = scalePolygon(x, 0.85f);
        // m_blocks[i] = scalePolygon(x, 0.85f);
    }
}
#pragma endregion

void CityGen::sweepToBlocks()
{
    cout<<"\n\n\n\n\nNEW\n\n";
    m_chunks = m_voronoiCells;
    m_blocks.clear();
    m_blocks.resize(m_chunks.size());

    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];
    
        m_chunks[i] = scalePolygon(x, 0.85f);
        // m_blocks[i] = scalePolygon(x, 0.85f);
    }

  // Cuts chunks into blocks
    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];

        pair<Point,Point> largestPointPair = findLargestEdge(x);

        Line largestEdge = CreateLineFromPoints(largestPointPair.first,largestPointPair.second);

        Point midpoint = findMidpoint(largestPointPair.first,largestPointPair.second);
        addDebug1(midpoint);
        
        Point centerPoint = getCentroid(x);
        addDebug1(centerPoint);
        
        Point directionVector = getDirectionVector(midpoint,centerPoint);
        addLineDebug({0,0},directionVector);

        Point goToPoint = midpoint;
        largestEdge = moveToPoint(largestEdge,goToPoint);

        int moveAmount = 10.0f;
        int iterationAmount = 0;
        float minEdge = 2.0f;
        bool fuck = false;

        // float pArea = calculatePolygonArea(clipped.first);
        // float nArea = calculatePolygonArea(clipped.second);

        double eval = evaluate(largestEdge, centerPoint);
        bool keepPositiveSide = eval > 0;

        cout<<"b4 - keep Positive Side: "<<keepPositiveSide<<endl;

        while(true)
        {
            cout <<"\n"<<iterationAmount<< " -Before move: Line: " << largestEdge.a << "x + " << largestEdge.b << "y + " << largestEdge.c << " = 0" << endl;
     
            cout << "goToPoint before: (" << goToPoint.x << ", " << goToPoint.y << ")\n";
            goToPoint = movePointInDirection(goToPoint, directionVector, moveAmount);
            cout << "goToPoint after: (" << goToPoint.x << ", " << goToPoint.y << ")\n";
            addDebug(goToPoint);

            largestEdge = moveToPoint(largestEdge,goToPoint);
            cout << "After move: Line: " << largestEdge.a << "x + " << largestEdge.b << "y + " << largestEdge.c << " = 0" << endl;

            debugs.push_back(findIntersectionsWithBoundary(largestEdge));

            iterationAmount++;
            auto clipped = clipPolygon(x, largestEdge);
            vector<Point> positive = clipped.first;
            vector<Point> negative = clipped.second;

            if (positive.empty() || negative.empty()){
                cout << "Clipping failed: Positive size: " << positive.size()
                        << ", Negative size: " << negative.size() << "\n";
                addDebug1(goToPoint); // Visualize the midpoint and sweep line
                debugs.push_back(findIntersectionsWithBoundary(largestEdge));
                debugs1.push_back(positive);
                debugs1.push_back(negative);
                fuck = true;
            }

            if((positive.empty() && negative.empty()))
                {cout<<"1 termination 1"<<endl; break;}
            
            // float drawPos = findSmallestEdgeAmount(positive);
            // float drawNeg = findSmallestEdgeAmount(negative);
            // cout<<"Smallest Positive Edge: "<<drawPos
            // <<"\nSmallest Negative Edge: "<<drawNeg<<endl;
            // if((drawPos <= minEdge && positive.size() <= 4) 
            // || ( drawNeg <= minEdge && negative.size() <= 4))

            bool drawPos = findSmallestEdgeAmount(positive, minEdge);
            bool drawNeg = findSmallestEdgeAmount(negative, minEdge);

            if((drawPos && positive.size() <= 4) 
            || ( drawNeg && negative.size() <= 4))
            {
                cout<<"---a"<<endl;
                m_blocks.push_back(x);
                // x.clear();
                break;
            }
            else if(keepPositiveSide)
            {
                cout<<"---b\n"<<endl;
                if(positive.empty()) {cout<<"2 termination 2"<<endl; break;}

                m_blocks.push_back(negative);
                x = positive;
            }else
            {
                cout<<"---c\n"<<endl;
                if(negative.empty()) {cout<<"3 termination 3"<<endl; break;}

                m_blocks.push_back(positive);
                x = negative;
            }
            if(fuck) return;
        }

        float distance = iterationAmount * moveAmount;
    }
}
void CityGen::buildVertexData() {
    /*
    for (size_t i = 0; i < m_blocks.size(); ++i)
    {
        const vector<Point>& cell = m_blocks[i];

        for (size_t j = 0; j < cell.size(); ++j) {

            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            CreateCubesAlongLine(p1,p2,2.0f);
        }
    }
    */
    // /*
    m_vertices.clear();

    vector<GLubyte> colors = {
    // 255, 0, 0,     // Bright Red
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

    // vector<GLubyte> colors = {
    //     255, 0, 0,   // Red
    //     0, 255, 0,   // Green
    //     0, 0, 255,   // Blue
    //     255, 255, 0, // Yellow
    //     255, 0, 255, // Magenta
    //     0, 255, 255, // Cyan
    //     128, 0, 128, // Purple
    //     255, 165, 0, // Orange
    //     0, 128, 0,   // Dark Green
    //     128, 128, 128 // Gray
    // };
    for (size_t i = 0; i < m_blocks.size(); ++i) {


        const vector<Point>& cell = m_blocks[i];
        
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        GLubyte a = 255;
        
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            Vertex v1 = { static_cast<GLfloat>(p1.x), 0.0f, static_cast<GLfloat>(p1.y), r, g, b, a };
            Vertex v2 = { static_cast<GLfloat>(p2.x), 0.0f, static_cast<GLfloat>(p2.y), r, g, b, a };

            m_vertices.push_back(v1);
            m_vertices.push_back(v2);
        }   
    }
    for (size_t j = 0; j < debugs.size(); ++j)
    {
        if(debugs[j].empty()) continue;
        Point p1 = debugs[j][0];
        Point p2 = debugs[j][1];
        Vertex v1 = { static_cast<GLfloat>(p1.x), 10.0f, static_cast<GLfloat>(p1.y), 255, 0, 0, 255};
        Vertex v2 = { static_cast<GLfloat>(p2.x), 10.0f, static_cast<GLfloat>(p2.y), 255, 0, 0, 255};
        m_vertices.push_back(v1);
        m_vertices.push_back(v2);
    }
    for (size_t j = 0; j < debugs1.size(); ++j)
    {
        if(debugs1[j].empty()) continue;
        Point p1 = debugs1[j][0];
        Point p2 = debugs1[j][1];
        Vertex v1 = { static_cast<GLfloat>(p1.x), 20.0f, static_cast<GLfloat>(p1.y), 0, 0, 255, 255};
        Vertex v2 = { static_cast<GLfloat>(p2.x), 20.0f, static_cast<GLfloat>(p2.y), 0, 0, 255, 255};
        m_vertices.push_back(v1);
        m_vertices.push_back(v2);
    }
    for (size_t j = 0; j < debugs2.size(); ++j)
    {
        if(debugs2[j].empty()) continue;
        Point p1 = debugs2[j][0];
        Point p2 = debugs2[j][1];
        Vertex v1 = { static_cast<GLfloat>(p1.x), 30.0f, static_cast<GLfloat>(p1.y), 0, 255, 0, 255};
        Vertex v2 = { static_cast<GLfloat>(p2.x), 30.0f, static_cast<GLfloat>(p2.y), 0, 255, 0, 255};
        m_vertices.push_back(v1);
        m_vertices.push_back(v2);
    }

    m_numVertices = static_cast<int>(m_vertices.size());// * 5.0f); // sussyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

    // Create and bind VAO and VBO
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size() * sizeof(Vertex), m_vertices.data(), GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(0);

    // Color attribute
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (void*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);

    // Unbind VAO
    glBindVertexArray(0);
    // */
}

void CityGen::generate(GLuint program) {
    m_program = program;
    
    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    m_sites = generateSites(numSites, seed);

    // m_sites = {
    // {-50,-50},
    // {50,50}};

    computeVoronoiDiagram();
    computeChunks();
    sweepToBlocks();
    buildVertexData();
}

#pragma region Render
void CityGen::deGenerate() {
    // Clean up OpenGL resources
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
    cubes.clear();
    m_voronoiCells.clear();
}
void CityGen::reGenerate() {
    deGenerate();
    generate(m_program);
}
void CityGen::render(const glm::mat4& proj, const glm::mat4& view, float m_timer) {
    glUseProgram(m_program);

    // Set uniform variables
    GLint projLoc = glGetUniformLocation(m_program, "projection");
    GLint viewLoc = glGetUniformLocation(m_program, "view");
    GLint worldLoc = glGetUniformLocation(m_program, "world");
    GLint timeLoc = glGetUniformLocation(m_program, "u_time");

    if (projLoc != -1) glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));
    if (viewLoc != -1) glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    if (worldLoc != -1) glUniformMatrix4fv(worldLoc, 1, GL_FALSE, glm::value_ptr(glm::mat4(1.0f)));
    if (timeLoc != -1) glUniformMatrix4fv(timeLoc, 1, GL_FALSE, glm::value_ptr(glm::mat4(m_timer)));

    glLineWidth(3.0f); // Set the line thickness to 3.0 (default is 1.0)
    glEnable(GL_LINE_SMOOTH);
    // glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    // Bind VAO and draw lines
    glBindVertexArray(m_vao);
    glDrawArrays(GL_LINES, 0, m_numVertices);
    glBindVertexArray(0);
    
    glUseProgram(0);
}
#pragma endregion