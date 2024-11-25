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

    double slope = dy / dx;

    double a = -slope;
    double b = 1.0;
    double c = -(a * p1.x + b * p2.y);

    return {a, b, c};

}
Line CityGen::perpendicularBisector(const Point& p1, const Point& p2) {
    // Midpoint
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;

    // Slope of the line between p1 and p2
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    if (std::abs(dx) < 1e-9) {
        // Vertical line; perpendicular bisector is horizontal
        return {0.0, 1.0, -my};
    }

    double slope = dy / dx;
    // Slope of the perpendicular bisector
    double perpSlope = -dx / dy;

    // Line equation: (y - my) = perpSlope * (x - mx)
    // Rearranged to a*x + b*y + c = 0
    double a = -perpSlope;
    double b = 1.0;
    double c = -(a * mx + b * my);

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
#pragma endregion
#pragma region Debug
void CityGen::CreateCubesAlongLine(const Point& start, const Point& end, int height) {
    glm::vec3 direction(end.x - start.x, 0.0f, end.y - start.y);

    glm::vec3 normalizedDir = glm::normalize(direction);

    float lineLength = glm::length(direction);
    float spacing = 5.0f;

    for (int i = 0; i < lineLength/spacing; ++i) {
        glm::vec3 position = glm::vec3(start.x, 0.0f, start.y) + normalizedDir * (spacing * i);
        
        Cube c(position, vec3(3.0f,3.0f,height), vec3(0),vec4(1));
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

Line moveLinePerpendicularly(const Line& line, double distance) {
    double magnitude = sqrt(line.a * line.a + line.b * line.b);
    double newC = line.c + distance * magnitude; // Shift the line by the given distance
    return {line.a, line.b, newC};
}
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


void CityGen::computeChunks() {
    m_voronoiCells;
    m_chunks.clear();
    m_chunks.resize(m_voronoiCells.size());
    m_blocks.clear();
    m_blocks.resize(m_chunks.size());

    for (size_t i = 0; i < m_voronoiCells.size(); i++) {
        vector<Point> chunk = m_voronoiCells[i];
        m_chunks.push_back(chunk);

        // Line cut = CreateLineFromPoints(chunk[0],chunk[1]);
        // cut = moveLinePerpendicularly(cut, 20.0f);
        // auto clipped = clipPolygon(chunk, cut);
        // vector<Point> positive = clipped.first;
        // vector<Point> negative = clipped.second;
        // if (!positive.empty()) {
        //     m_chunks.push_back(positive);
        // }
        // if (!negative.empty()) {
        //     m_chunks.push_back(negative);
        // }
        m_blocks[i] = scalePolygon(chunk, 0.85f);
        // m_blocks[i] = offsetPolygonInward(chunk, 10.0f);
    }

    // for (size_t i = 0; i < m_chunks.size(); i++)
    // {
    //     vector<Point> chunk = offsetPolygonInward(m_chunks[i], 5.0f);
    //     m_blocks[i] = chunk;
    //     m_blocks[i] = offsetPolygonInward(chunk, 5.0f);
    // }

}
void CityGen::sweepToBlocks()
{
    m_voronoiCells;
    m_chunks;
    m_blocks.clear();

    for (size_t i = 0; i < m_chunks.size(); i++) {

        vector<Point> chunk = m_chunks[i];

        Line cut = CreateLineFromPoints(chunk[0],chunk[1]);
        
        cut = moveLinePerpendicularly(cut, 40.0f);

        auto clipped = clipPolygon(chunk, cut);

        vector<Point> positive = clipped.first;
        vector<Point> negative = clipped.second;

        if (!positive.empty()) {
            m_blocks.push_back(positive);
        }
        if (!negative.empty()) {
            m_blocks.push_back(negative);
        }
    }
}
void CityGen::buildVertexData() {
    for (size_t i = 0; i < m_blocks.size(); ++i)
    {
        const vector<Point>& cell = m_blocks[i];

        for (size_t j = 0; j < cell.size(); ++j) {

            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            CreateCubesAlongLine(p1,p2,2.0f);
        }
    }
    /*m_vertices.clear();

    vector<GLubyte> colors = {
        255, 0, 0,   // Red
        0, 255, 0,   // Green
        0, 0, 255,   // Blue
        255, 255, 0, // Yellow
        255, 0, 255, // Magenta
        0, 255, 255, // Cyan
        128, 0, 128, // Purple
        255, 165, 0, // Orange
        0, 128, 0,   // Dark Green
        128, 128, 128 // Gray
    };
    for (size_t i = 0; i < m_blocks.size(); ++i) {

        const vector<Point>& cell = m_blocks[i];
        
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        GLubyte a = 255;
        
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            Vertex v1 = { static_cast<GLfloat>(p1.x), static_cast<GLfloat>(p1.y), 0.0f, r, g, b, a };
            Vertex v2 = { static_cast<GLfloat>(p2.x), static_cast<GLfloat>(p2.y), 0.0f, r, g, b, a };

            m_vertices.push_back(v1);
            m_vertices.push_back(v2);
        }
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
    */
}

void CityGen::generate(GLuint program) {
    m_program = program;
    
    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    m_sites = generateSites(numSites, seed);

    // for (const auto& site : m_sites) {
    //     Cube c(vec3(site.x,0.0f,site.y), vec3(0.0f), vec3(0), vec4(1.0f));
    //     c.init(m_program);
    //     cubes.push_back(c);
    // }

    computeVoronoiDiagram();
    computeChunks();
    // sweepToBlocks();
    buildVertexData();

}

#pragma Render
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
    GLint modelLoc = glGetUniformLocation(m_program, "model");

    glm::mat4 model = glm::mat4(1.0f);
    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    // Bind VAO and draw lines
    glBindVertexArray(m_vao);
    // glDrawArrays(GL_LINES, 0, m_numVertices);
    glBindVertexArray(0);

    for (auto& c : cubes)
    {
        c.draw(proj, view, m_timer);
    }
    
    glUseProgram(0);
}
#pragma endregion