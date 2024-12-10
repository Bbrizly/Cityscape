#include "PolygonUtils.h"

double PolygonUtils::distanceBetweenPoints(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
}

std::pair<Point,Point> PolygonUtils::findLargestEdge(std::vector<Point> polygon)
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
        const Point& p2 = polygon[(i + 1) % n];
        double length = distanceBetweenPoints(p1,p2);
        if (length > maxLength) {
            maxLength = length;
            largestEdge = {p1, p2};
        }
    }
    return largestEdge;
}

std::pair<std::vector<Point>, std::vector<Point>> PolygonUtils::clipPolygon(std::vector<Point> polygon, Line l) {
    std::vector<Point> polygonPositive;
    std::vector<Point> polygonNegative;
    if (polygon.empty()) {
        return {polygonPositive, polygonNegative};
    }
    Point prev = polygon.back();
    double prevEval = GeometryUtils::evaluate(l, prev);
    for (const Point& curr : polygon) {
        double currEval = GeometryUtils::evaluate(l, curr);
        bool prevInside = prevEval >= 0;
        bool currInside = currEval >= 0;
        if (currInside) {
            if (!prevInside) {
                Point intersection;
                if (GeometryUtils::lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    polygonPositive.push_back(intersection);
                    polygonNegative.push_back(intersection);
                }
            }
            polygonPositive.push_back(curr);
        }
        else {
            if (prevInside) {
                Point intersection;
                if (GeometryUtils::lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    polygonPositive.push_back(intersection);
                    polygonNegative.push_back(intersection);
                }
            }
            polygonNegative.push_back(curr);
        }
        prev = curr; 
        prevEval = currEval;
    }
    return {polygonPositive, polygonNegative};
}

Point PolygonUtils::getCentroid(std::vector<Point> polygon) {
    if (polygon.empty()) return {};
    double centroidX = 0.0;
    double centroidY = 0.0;
    for (const Point& p : polygon) {
        centroidX += p.x + 0.5f;
        centroidY += p.y + 0.5f;
    }
    centroidX /= polygon.size();
    centroidY /= polygon.size();
    return {centroidX, centroidY};
}

Line PolygonUtils::moveLineToCenter(Line l, std::vector<Point> polygon) {
    Point centroid = getCentroid(polygon);
    return GeometryUtils::moveToPoint(l, centroid);
}

std::vector<Point> PolygonUtils::scalePolygon(const std::vector<Point>& polygon, float scaleFactor) {
    if (polygon.empty()) return {};
    double centroidX = 0.0;
    double centroidY = 0.0;
    for (const Point& p : polygon) {
        centroidX += p.x;
        centroidY += p.y;
    }
    centroidX /= polygon.size();
    centroidY /= polygon.size();

    std::vector<Point> scaledPolygon;
    for (const Point& p : polygon) {
        double newX = centroidX + (p.x - centroidX) * scaleFactor;
        double newY = centroidY + (p.y - centroidY) * scaleFactor;
        scaledPolygon.push_back({newX, newY});
    }
    return scaledPolygon;
}

float PolygonUtils::calculatePolygonArea(const std::vector<Point>& polygon) {
    float area = 0.0f;
    size_t n = polygon.size();
    for (size_t i = 0; i < n; i++) {
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % n];
        area += (p1.x * p2.y - p2.x * p1.y);
    }
    return std::abs(area) / 2.0f;
}

bool PolygonUtils::findSmallestEdgeAmount(std::vector<Point> polygon, float minEdgeLength) {
    if (polygon.size() < 2) {
        return false;
    }
    size_t n = polygon.size();
    int edgesBelowMin = 0;
    for (size_t i = 0; i < n; ++i) {
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[(i + 1) % n];
        double length = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
        if (length < minEdgeLength) {
            edgesBelowMin++;
            if (edgesBelowMin >= 2) {
                return true;
            }
        }
    }
    return false;
}

std::vector<Vertex> PolygonUtils::fanTriangulatePolygon(const std::vector<Point>& polygon, 
                                                   const glm::vec3& normal, float height, float layer,
                                                   GLubyte r, GLubyte g, GLubyte b, GLubyte a)
{
    std::vector<Vertex> result;
    if (polygon.size() < 3) return result;

    double minX = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    double minZ = std::numeric_limits<double>::max();
    double maxZ = -std::numeric_limits<double>::max();

    for (const auto& p : polygon) {
        if (p.x < minX) minX = p.x;
        if (p.x > maxX) maxX = p.x;
        if (p.y < minZ) minZ = p.y;
        if (p.y > maxZ) maxZ = p.y;
    }
    float div = 4.0f;
    if(layer == 7) //if layer is sidewalk, add more divisions for scalability
    {
        div = 10.0f;
    }

    double dx = (maxX - minX) / div;
    double dz = (maxZ - minZ) / div;
    if (dx < 1e-9) dx = 1.0; 
    if (dz < 1e-9) dz = 1.0;


    std::random_device rd;
    std::mt19937 gen(rd());
    // std::uniform_real_distribution<float> scaleDist(uvMinScale, uvMaxScale);
    std::uniform_real_distribution<float> rotDist(0.0f, 360.0f);

    float angle = glm::radians(rotDist(gen));

    float ca = cos(angle);
    float sa = sin(angle);

    std::vector<Vertex> tempVerts;
    tempVerts.reserve(polygon.size());

    for (const auto& p : polygon) {
        Vertex v;
        v.x = (GLfloat)p.x;
        v.y = (GLfloat)height;
        v.z = (GLfloat)p.y;
        v.r = r; v.g = g; v.b = b; v.a = a;
        v.nx = normal.x; v.ny = normal.y; v.nz = normal.z;
        v.layer = layer;

        // UV based on bounding box:
        float u = (GLfloat)((p.x - minX)/dx);
        float vcoord = (GLfloat)((p.y - minZ)/dz);

        // Rotate
        float uRot = u * ca - vcoord * sa;
        float vRot = u * sa + vcoord * ca;

        v.u = uRot;
        v.v = vRot;

        tempVerts.push_back(v);
    }

    for (size_t i = 1; i < polygon.size() - 1; i++) {
        result.push_back(tempVerts[0]);
        result.push_back(tempVerts[i+1]);
        result.push_back(tempVerts[i]);
    }

    return result;
}
