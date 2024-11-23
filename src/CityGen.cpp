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

vector<Building> buildings;

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

void CityGen::CreateBuildingsAlongLine(const Point& start, const Point& end, int height) {
    // Calculate the direction vector between start and end points
    glm::vec3 direction(end.x - start.x, 0.0f, end.y - start.y);

    // Normalize the direction vector
    glm::vec3 normalizedDir = glm::normalize(direction);

    // Calculate spacing between buildings
    float lineLength = glm::length(direction);
    float spacing = 5.0f;
    // float spacing = lineLength / static_cast<float>(numBuildings - 1);

    for (int i = 0; i < lineLength/spacing; ++i) {
        glm::vec3 position = glm::vec3(start.x, 0.0f, start.y) + normalizedDir * (spacing * i);

        Building b(position, height);
        b.init(m_program);
        buildings.push_back(b);

    }
}

vector<Point> CityGen::findIntersectionsWithBoundary(const Line& line) {
    vector<Point> intersections;

    // Define boundary lines
    Line top = {0, 1, -maxY};          // Horizontal: y = maxY
    Line bottom = {0, 1, -minY};       // Horizontal: y = minY
    Line left = {1, 0, -minX};         // Vertical: x = minX
    Line right = {1, 0, -maxX};        // Vertical: x = maxX

    // Check intersection with each boundary line
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
    Building b(vec3(intersections[1].x,0.0f,intersections[1].y),5.0f);
    b.init(m_program);
    buildings.push_back(b);
    
    Building b1(vec3(intersections[0].x,0.0f,intersections[0].y),5.0f);
    b1.init(m_program);
    buildings.push_back(b1);

    return intersections;
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

double CityGen::evaluate(const Line& l, const Point& p) {
    return (l.a * p.x + l.b * p.y + l.c);
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

vector<Point> CityGen::clipPolygon(const vector<Point>& polygon, const Line& l, bool keepPositiveSide) 
{
    vector<Point> newPolygon;
    if (polygon.empty()) {
        cout << "\nPolygon is empty.";
        return newPolygon;
    }

    Point prev = polygon.back();
    double prevEval = evaluate(l, prev);

    for (const Point& curr : polygon) {
        double currEval = evaluate(l, curr);
        if (keepPositiveSide) {
            if (currEval >= 0) {
                if (prevEval < 0) {
                    //from outside to inside
                    Point intersection;
                    if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                        newPolygon.push_back(intersection);
                    }
                }
                //point is inside
                newPolygon.push_back(curr);
            } else if (prevEval >= 0) {
                //from inside to outside
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    newPolygon.push_back(intersection);
                }
            }
        } else {
            // If keeping the negative side
            if (currEval <= 0) {
                if (prevEval > 0) {
                    //from outside to inside
                    Point intersection;
                    if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                        newPolygon.push_back(intersection);
                    }
                }
                //point is inside
                newPolygon.push_back(curr);
            } else if (prevEval <= 0) {
                //inside to outside
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    newPolygon.push_back(intersection);
                }
            }
        }

        /*
        if (currEval >= 0) {
            if (prevEval < 0) {
                // Edge crosses the line from outside to inside
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    newPolygon.push_back(intersection);
                }
            }
            // Current point is inside
            newPolygon.push_back(curr);
        } else if (prevEval >= 0) {
            // Edge crosses the line from inside to outside
            Point intersection;
            if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                newPolygon.push_back(intersection);
            }
        }
        */
        prev = curr;
        prevEval = currEval;
    }
    return newPolygon;
    /*
    vector<Point> newPolygon;
    if (polygon.empty())
    {cout<<"\nPlygon is empty WTF?"; return newPolygon;}

    Point prev = polygon.back();
    double prevEval = evaluate(l, prev);

    for (const Point& curr : polygon) {
        double currEval = evaluate(l, curr);

        cout<<"\nPoint: x:"<<curr.x<<" y:"<<curr.y<<" Eval:"<<currEval<<endl;

        if (currEval >= 0) {
            if (prevEval < 0) {
                // Edge crosses the line from outside to inside
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    newPolygon.push_back(intersection);
                    CreateBuildingsAlongLine(intersection,prev,10.0f);  
                }
            }
            // Current point is inside
            newPolygon.push_back(curr);
        } else if (prevEval >= 0) {
            // Edge crosses the line from inside to outside
            Point intersection;
            if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                newPolygon.push_back(intersection);
                CreateBuildingsAlongLine(prev,intersection,10.0f);  

            }
        }
        prev = curr;
        prevEval = currEval;
    }
    
    for (size_t j = 0; j < polygon.size(); ++j) {
        Point p1 = polygon[j];
        Point p2 = polygon[(j + 1) % polygon.size()];

        CreateBuildingsAlongLine(p1,p2,10.0f);
    }
    
    return newPolygon;
    */
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
            
            cell = clipPolygon(cell, bisector, keepPositiveSide);

            if (cell.empty()) {
                cout<<"\nCell is empty somefuckinghow";
                break;
            }
        }
        cout<<"\n\n\n"<<endl;
        m_voronoiCells[i] = cell;
    }
}

void CityGen::buildVertexData() {
    m_vertices.clear();

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

    for (size_t i = 0; i < m_voronoiCells.size(); ++i) {

        const vector<Point>& cell = m_voronoiCells[i];
        
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        GLubyte a = 255;
        
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            CreateBuildingsAlongLine(p1,p2,10.0f);

            Vertex v1 = { static_cast<GLfloat>(p1.x), static_cast<GLfloat>(p1.y), 0.0f, r, g, b, a };
            Vertex v2 = { static_cast<GLfloat>(p2.x), static_cast<GLfloat>(p2.y), 0.0f, r, g, b, a };

            // m_vertices.push_back(v1);
            // m_vertices.push_back(v2);
        }
    }

    m_numVertices = static_cast<int>(m_vertices.size() * 3.0f); // sussyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

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
}

void CityGen::generate(GLuint program) {
    m_program = program;

    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    m_sites = generateSites(numSites, seed);

    cout << "Generated Sites:" << endl;
    for (const auto& site : m_sites) {
        cout << "Site: x = " << site.x << ", y = " << site.y << endl;
        Building b(vec3(site.x, 0.0f, site.y), 1000.0f);
        b.init(m_program);
        buildings.push_back(b);
    }

    computeVoronoiDiagram();
    buildVertexData();
}

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
    buildings.clear();
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
    glDrawArrays(GL_LINES, 0, m_numVertices);
    glBindVertexArray(0);

    for (auto& building : buildings)
    {
        building.draw(proj, view, m_timer);
    }
    
    glUseProgram(0);
}