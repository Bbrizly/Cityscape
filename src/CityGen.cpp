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

vector<Point> CityGen::generateSites(int numSites, unsigned int seed) {

    mt19937 gen(seed);
    uniform_real_distribution<double> distX(minX + 0, maxX - 0);
    uniform_real_distribution<double> distY(minY + 0, maxY - 0);

    vector<Point> sites;
    for (int i = 0; i < numSites; ++i) {
        float x = distX(gen);
        float y = distY(gen);
        sites.push_back({x, y});
    }
    return sites;
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

    // Avoid division by zero for vertical lines
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
    /*
    // Midpoint
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;

    {mx,my;}

    //slope of p1 to p2
    //y2-y1/x2-x1
    double slope = (p2.y - p1.y) / (p2.x - p1.x);

    //reciprocal = - 1/n
    double reciprocalSlope = -(1.0 / slope);

    // y - midpointy = slope * (x- midpointx)

    double a = -reciprocalSlope;
    double b = 1;
    double c = reciprocalSlope * mx - my;
        
    Building b0(vec3(mx, 0.0f, my),1.0f);
    b0.init(m_program);
    buildings.push_back(b0);

    vector<Point> intersections = findIntersectionsWithBoundary({a,b,c});    

    // CreateBuildingsAlongLine(intersections[0], intersections[1], 2);

    return {a,b,c};
    */
}

double CityGen::evaluate(const Line& l, const Point& p) {
    
    
    return (l.a * p.x + l.b * p.y + l.c);
}

bool CityGen::lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l1, Point& intersection) {
    p1; p2; l1;
    //create line 2 from p1 and p2

    double slope = (p2.y - p1.y) / (p2.x - p1.x);
    Line l2;
    l2.a = -slope;
    l2.b = 1;
    l2.c = slope * p1.x - p1.y;
    bool ret = findLineIntersection(l1, l2, intersection);

    return (ret);

    // Line l2;
    // l2.a = p2.y - p1.y; //dy
    // l2.b = p1.x - p2.x; //-dx
    // l2.c = -(l2.a * p1.x + l2.b * p1.y);
    /*double determinant = l1.a * l2.b - l2.a * l1.b;
    if (determinant == 0) {
        std::cout << "The lines are parallel or coincident, no unique intersection." << std::endl;
        return false;
    }

    // Calculate the intersection point
    intersection.x = (l1.b * l2.c - l2.b * l1.c) / determinant;
    intersection.y = (l2.a * l1.c - l1.a * l2.c) / determinant;

    
    Building b0(vec3(intersection.x, 0.0f, intersection.y),25);
    b0.init(m_program);
    buildings.push_back(b0);

    // if(fabs(determinant) < numeric_limits<double>::epsilon())
    //     {return false;}
    
    // intersection.x = (l2.b * -l.c - l.b * -l2.c) / determinant;
    // intersection.y = (l.a * -l2.c - l2.a * -l.c) / determinant;

    return true;
    */
    /*double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    double denominator = l.a * dx + l.b * dy;
    if (std::abs(denominator) < epsilon) {
        // Line segment is parallel to the line
        return false;
    }

    // double t = - (l.a * p1.x + l.b * p1.y + l.c) / denominator;
    double t = - evaluate(l,p1) / denominator;
    if (t < -epsilon || t > 1.0 + epsilon) {
        // Intersection is outside the segment
        return false;
    }

    intersection.x = p1.x + t * dx;
    intersection.y = p1.y + t * dy;

    return true;
    */
}

vector<Point> CityGen::clipPolygon(const vector<Point>& polygon, const Line& l) 
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
    // Initialize Voronoi cells
    m_voronoiCells.clear();
    m_voronoiCells.resize(m_sites.size());

    // Bounding rectangle as initial polygon
    vector<Point> boundingPolygon = {
        { minX, minY },
        { maxX, minY },
        { maxX, maxY },
        { minX, maxY }
    };

    // Compute Voronoi cell for each site
    for (size_t i = 0; i < m_sites.size(); i++) {
        if(i>0) return;
        const Point& site = m_sites[i];
        vector<Point> cell = boundingPolygon;

        for (size_t j = 0; j < m_sites.size(); j++) {
            if (i == j) continue;

            const Point& otherSite = m_sites[j];

            Line bisector = perpendicularBisector(site, otherSite);

            // Clip the cell polygon with the half-plane
            cout<<"Cell: "<<i<<", How many did it survive?: "<<j<<endl;
            
            // cell = clipPolygon(cell, bisector);
            /* Testing removing the clip polygon code*/
            vector<Point> newPolygon;
            vector<Point> polygon = cell;
            Line l = bisector;

            if (polygon.empty()) {
                cout << "\nPolygon is empty.";
                continue;}

            Point prev = polygon.back();
            double prevEval = evaluate(l, prev);

            for (const Point& curr : polygon) {
                double currEval = evaluate(l, curr);

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
                prev = curr;
                prevEval = currEval;
            }

            if (cell.empty()) {
                cout<<"\nCell is empty somefuckinghow";
                break; // The cell is empty, no need to continue
            }
        }
        cout<<"\n\n\n"<<endl;
        m_voronoiCells[i] = cell;
    }
}

void CityGen::buildVertexData() {
    m_vertices.clear();

    // Assign colors for different cells
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

    // Build vertices for rendering edges
    for (size_t i = 0; i < m_voronoiCells.size(); ++i) {

        const vector<Point>& cell = m_voronoiCells[i];
        
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        GLubyte a = 255;
        
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            // CreateBuildingsAlongLine(p1,p2,10.0f);

            Vertex v1 = { static_cast<GLfloat>(p1.x), static_cast<GLfloat>(p1.y), 0.0f, r, g, b, a };
            Vertex v2 = { static_cast<GLfloat>(p2.x), static_cast<GLfloat>(p2.y), 0.0f, r, g, b, a };

            // Vertex v1 = { static_cast<GLfloat>(p1.x), 0.0f, static_cast<GLfloat>(p1.y), r, g, b, a };
            // Vertex v2 = { static_cast<GLfloat>(p2.x), 0.0f, static_cast<GLfloat>(p2.y), r, g, b, a };

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

    Building b(vec3(0.0f),0.01f);
    b.init(m_program);
    buildings.push_back(b);
    

    // Generate random sites
    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    // m_sites = generateSites(numSites, seed);

    m_sites = {
        {0.0, 0.0},
        {50.0, 50.0},
        {-50.0, 50.0}//,
        // {-50.0, -50.0},
        // {50.0, 50.0},
        // {25.0, 25.0}
    };

    cout << "Generated Sites:" << endl;

    for (const auto& site : m_sites) {
        cout << "Site: x = " << site.x << ", y = " << site.y << endl;
        Building b(vec3(site.x, 0.0f, site.y), 1000.0f);
        b.init(m_program);
        buildings.push_back(b);
    }

    // Compute Voronoi diagram
    computeVoronoiDiagram();

    // Build vertex data for rendering
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

/*
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
    uniform_real_distribution<double> distX(minX + 0, maxX - 0);
    uniform_real_distribution<double> distY(minY + 0, maxY - 0);

    vector<Point> sites;
    for (int i = 0; i < numSites; ++i) {
        float x = distX(gen);
        float y = distY(gen);
        sites.push_back({x, y});
        
        Building b(vec3(x, 0.0f, y),1000.0f);
        b.init(m_program);
        buildings.push_back(b);

    }
    return sites;
}

Line CityGen::perpendicularBisector(const Point& p1, const Point& p2) {
    // Midpoint
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;

    // Direction vector (p1 to p2)
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    // Perpendicular direction
    double a = -dy;
    double b = dx;

    // Line equation: a(x - mx) + b(y - my) = 0 => ax + by + c = 0
    double c = - (a * mx + b * my);

    return { a, b, c };
}

double CityGen::evaluate(const Line& l, const Point& p) {
    return l.a * p.x + l.b * p.y + l.c;
}

bool CityGen::lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    double denominator = l.a * dx + l.b * dy;
    if (std::abs(denominator) < epsilon) {
        // Line segment is parallel to the line
        return false;
    }

    // double t = - (l.a * p1.x + l.b * p1.y + l.c) / denominator;
    double t = - evaluate(l,p1) / denominator;
    if (t < -epsilon || t > 1.0 + epsilon) {
        // Intersection is outside the segment
        return false;
    }

    intersection.x = p1.x + t * dx;
    intersection.y = p1.y + t * dy;

    // // Building b(vec3(p1.x, 0.0f, p1.y),0.25f);
    // Building b(vec3(intersection.x, 0.0f, intersection.y),0.1f);
    // b.init(cubeShader);
    // buildings.push_back(b);
    return true;
}

vector<Point> CityGen::clipPolygon(const vector<Point>& polygon, const Line& l, const Point& site) {
    vector<Point> newPolygon;
    if (polygon.empty()) return newPolygon;

    Point prev = polygon.back();
    double prevEval = evaluate(l, prev);
    double prevSide = prevEval >= -epsilon ? 1 : -1;

    for (const Point& curr : polygon) {
        double currEval = evaluate(l, curr);
        double currSide = currEval >= -epsilon ? 1 : -1;

        if (currSide >= 0) {
            if (prevSide < 0) {
                // Edge crosses the line from outside to inside
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    newPolygon.push_back(intersection);
                }
            }
            // Current point is inside
            newPolygon.push_back(curr);
        } else if (prevSide >= 0) {
            // Edge crosses the line from inside to outside
            Point intersection;
            if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                newPolygon.push_back(intersection);
            }
        }
        prev = curr;
        prevSide = currSide;
    }
    return newPolygon;
}

void CityGen::computeVoronoiDiagram() {
    // Initialize Voronoi cells
    m_voronoiCells.clear();
    m_voronoiCells.resize(m_sites.size());

    // Bounding rectangle as initial polygon
    vector<Point> boundingPolygon = {
        { minX, minY },
        { maxX, minY },
        { maxX, maxY },
        { minX, maxY }
    };

    // Compute Voronoi cell for each site
    for (size_t i = 0; i < m_sites.size(); ++i) {
        const Point& site = m_sites[i];
        vector<Point> cell = boundingPolygon;
        
        // Building b(vec3(site.x, 0.0f, site.y),10.0f);
        // b.init(m_program);
        // buildings.push_back(b);

        for (size_t j = 0; j < m_sites.size(); ++j) {
            if (i == j) continue;

            const Point& otherSite = m_sites[j];

            Line bisector = perpendicularBisector(site, otherSite);

            // Clip the cell polygon with the half-plane
            cell = clipPolygon(cell, bisector, site);

            if (cell.empty()) {
                break; // The cell is empty, no need to continue
            }
        }
        m_voronoiCells[i] = cell;
        

        for (auto& point : cell)
        {
            // Building b(vec3(point.x, 0.0f, point.y),1000.0f);
            // b.init(m_program);
            // buildings.push_back(b);
        }
    }
}

void CityGen::CreateBuildingsAlongLine(const Point& start, const Point& end, int numBuildings) {
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

        Building b(position, 0.1f);
        b.init(m_program);
        buildings.push_back(b);

    }
}

void CityGen::buildVertexData() {
    m_vertices.clear();

    // Assign colors for different cells
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

    // Build vertices for rendering edges
    for (size_t i = 0; i < m_voronoiCells.size(); ++i) {

        const vector<Point>& cell = m_voronoiCells[i];
        
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        GLubyte a = 255;

        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];

            CreateBuildingsAlongLine(p1,p2,5);

            Vertex v1 = { static_cast<GLfloat>(p1.x), static_cast<GLfloat>(p1.y), 0.0f, r, g, b, a };
            Vertex v2 = { static_cast<GLfloat>(p2.x), static_cast<GLfloat>(p2.y), 0.0f, r, g, b, a };

            // Vertex v1 = { static_cast<GLfloat>(p1.x), 0.0f, static_cast<GLfloat>(p1.y), r, g, b, a };
            // Vertex v2 = { static_cast<GLfloat>(p2.x), 0.0f, static_cast<GLfloat>(p2.y), r, g, b, a };

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

    Building b1(vec3(0.0f,0.0f,10.0f),0.01f);
    b1.init(m_program);
    buildings.push_back(b1);
    
    Building b2(vec3(20.0f,0.0f,0.0f),0.01f);
    b2.init(m_program);
    buildings.push_back(b2);

    Building b(vec3(0.0f),0.01f);
    b.init(m_program);
    buildings.push_back(b);
    

    // Generate random sites
    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    m_sites = generateSites(numSites, seed);

    cout << "Generated Sites:" << endl;

    for (const auto& site : m_sites) {
        cout << "Site: x = " << site.x << ", y = " << site.y << endl;
        // Building b(vec3(site.x, 0.0f, site.y),1000.0f);
        // b.init(m_program);
        // buildings.push_back(b);
    }

    // Compute Voronoi diagram
    computeVoronoiDiagram();

    // Build vertex data for rendering
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
*/

/*
#include "CityGen.h"
#include "Building.h"
#include "Road.h"
#include "Cube.h"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <limits>
using namespace std;
using namespace glm;

vector<Building> buildings;
vector<Road> roads;

const int gridSize = 20;
const float spacing = 30.0f;
const float citySpacing = 1.0f;
const int sectorSize = 25;
const int skyScraperSector = 5;

const int numSites = 2;

struct Site {
    float x, z;
};

struct Point2D {
    float x, y;
};

struct Line {
    float a, b, c; // Line equation: ax + by + c = 0
};

struct Edge {
    Point2D p1, p2;
};

vector<Site> generateVoronoiSites(mt19937& gen)
{
    uniform_real_distribution<float> siteDist(-(gridSize * citySpacing), gridSize * citySpacing);
    vector<Site> sites;
    for (int i = 0; i < numSites; ++i) {
        sites.push_back({ siteDist(gen), siteDist(gen) });
    }
    return sites;
}

// Function to compute the perpendicular bisector between two points
Line computePerpendicularBisector(const Point2D& p1, const Point2D& p2) {
    // Midpoint
    Point2D mid;
    mid.x = (p1.x + p2.x) / 2.0f;
    mid.y = (p1.y + p2.y) / 2.0f;

    // Direction vector from p1 to p2
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;

    // Normal vector to the line (perpendicular to the direction vector)
    float a = -dy;
    float b = dx;

    // Line equation: a*x + b*y + c = 0
    float c = - (a * mid.x + b * mid.y);

    return {a, b, c};
}

// Function to evaluate the line equation at a point
float evaluateLine(const Line& l, const Point2D& p) {
    return l.a * p.x + l.b * p.y + l.c;
}

// Function to compute the intersection point between a line segment and a line
Point2D computeIntersection(const Point2D& S, const Point2D& E, const Line& l) {
    float denom = l.a * (E.x - S.x) + l.b * (E.y - S.y);
    const float epsilon = 1e-6f;
    if (std::abs(denom) < epsilon) {
        // Line and segment are parallel, return midpoint (could also handle differently)
        return {(S.x + E.x) / 2.0f, (S.y + E.y) / 2.0f};
    } else {
        float t = - (l.a * S.x + l.b * S.y + l.c) / denom;
        // No need to clamp t since we are intersecting with infinite lines
        Point2D I;
        I.x = S.x + t * (E.x - S.x);
        I.y = S.y + t * (E.y - S.y);
        return I;
    }
}

// Function to clip a convex polygon with a half-plane
std::vector<Point2D> clipPolygonWithHalfPlane(const std::vector<Point2D>& polygon, const Line& l, bool keepPositiveSide) {
    std::vector<Point2D> outputList;

    if (polygon.empty()) {
        return outputList;
    }

    Point2D S = polygon.back();
    float S_value = evaluateLine(l, S);
    for (const Point2D& E : polygon) {
        float E_value = evaluateLine(l, E);

        if ((E_value >= 0 && keepPositiveSide) || (E_value <= 0 && !keepPositiveSide)) {
            if ((S_value >= 0 && keepPositiveSide) || (S_value <= 0 && !keepPositiveSide)) {
                // Case 1: E and S are inside
                outputList.push_back(E);
            } else {
                // Case 4: E is inside, S is outside
                Point2D I = computeIntersection(S, E, l);
                outputList.push_back(I);
                outputList.push_back(E);
            }
        } else {
            if ((S_value >= 0 && keepPositiveSide) || (S_value <= 0 && !keepPositiveSide)) {
                // Case 2: E is outside, S is inside
                Point2D I = computeIntersection(S, E, l);
                outputList.push_back(I);
            } else {
                // Case 3: Both outside, do nothing
            }
        }
        S = E;
        S_value = E_value;
    }
    return outputList;
}

// Function to compute Voronoi cells
vector<vector<Point2D>> computeVoronoiCells(const vector<Site>& sites, float minX, float maxX, float minY, float maxY) {
    vector<vector<Point2D>> voronoiCells;
    // Define the bounding rectangle
    vector<Point2D> boundingPolygon = {
        {minX, minY},
        {maxX, minY},
        {maxX, maxY},
        {minX, maxY}
    };

    // Convert Site to Point2D
    vector<Point2D> points;
    for (const auto& site : sites) {
        points.push_back({site.x, site.z});
    }

    for (size_t i = 0; i < points.size(); ++i) {
        vector<Point2D> cell = boundingPolygon;
        Point2D site = points[i];

        for (size_t j = 0; j < points.size(); ++j) {
            if (i == j) continue;
            Point2D otherSite = points[j];
            Line bisector = computePerpendicularBisector(site, otherSite);

            // Determine which side to keep
            float siteValue = evaluateLine(bisector, site);
            float otherValue = evaluateLine(bisector, otherSite);
            bool keepPositiveSide = (siteValue > otherValue);

            // Clip the cell polygon with the half-plane
            cell = clipPolygonWithHalfPlane(cell, bisector, keepPositiveSide);

            if (cell.empty()) {
                break; // The cell is empty, no need to continue
            }
        }
        voronoiCells.push_back(cell);
    }
    return voronoiCells;
}

// Function to extract edges from Voronoi cells
vector<Edge> extractVoronoiEdges(const vector<vector<Point2D>>& voronoiCells) {
    vector<Edge> edges;
    for (const auto& cell : voronoiCells) {
        for (size_t i = 0; i < cell.size(); ++i) {
            Edge edge;
            edge.p1 = cell[i];
            edge.p2 = cell[(i + 1) % cell.size()];
            edges.push_back(edge);
        }
    }
    return edges;
}

// Function to create roads from Voronoi edges
void createRoadsFromVoronoiEdges(const vector<Edge>& edges, vector<Road>& roads, GLuint program) {
    for (const auto& edge : edges) {
        glm::vec3 midpoint = glm::vec3((edge.p1.x + edge.p2.x) / 2.0f, 0.0f, (edge.p1.y + edge.p2.y) / 2.0f);
        Road road(midpoint, 0);
        road.init(program);
        roads.push_back(road);

        Road road1(glm::vec3((edge.p1.x), 0.0f, (edge.p1.y)), 0);
        road1.init(program);
        roads.push_back(road1);

        Road road2(glm::vec3((edge.p2.x), 0.0f, (edge.p2.y)), 0);
        road2.init(program);
        roads.push_back(road2);
    }
}

// Function to place buildings within Voronoi cells
void placeBuildingsInVoronoiCells(const vector<vector<Point2D>>& voronoiCells, mt19937& gen, GLuint program, float scaleHeight) {
    uniform_real_distribution<float> heightDist(1.0f, scaleHeight);

    for (const auto& cell : voronoiCells) {
        // Compute the centroid of the cell
        float centroidX = 0.0f, centroidY = 0.0f;
        for (const auto& point : cell) {
            centroidX += point.x;
            centroidY += point.y;
        }
        centroidX /= cell.size();
        centroidY /= cell.size();

        // Place a building at the centroid
        vec3 buildingPosition = vec3(centroidX, 0.0f, centroidY);

        // Determine building height
        float height = heightDist(gen);

        Building b(buildingPosition, height);
        b.init(program);
        buildings.push_back(b);
    }
}

void CityGen::generate(GLuint program)
{
    m_program = program;

    random_device rd;
    mt19937 gen(rd());

    vector<Site> sites = generateVoronoiSites(gen);

    cout << "Voronoi Sites:" << endl;
    for (const auto& site : sites) {
        cout << "Site: x = " << site.x << ", z = " << site.z << endl;
        Building b(vec3(site.x, 0.0f, site.z), 1.0f);
        b.init(program);
        buildings.push_back(b);
    }

    // Define the bounding rectangle dimensions
    float minX = - (gridSize * citySpacing);
    float maxX = gridSize * citySpacing;
    float minY = - (gridSize * citySpacing);
    float maxY = gridSize * citySpacing;

    // Compute Voronoi cells
    vector<vector<Point2D>> voronoiCells = computeVoronoiCells(sites, minX, maxX, minY, maxY);

    // Extract edges from Voronoi cells
    vector<Edge> voronoiEdges = extractVoronoiEdges(voronoiCells);

    // Create roads from Voronoi edges
    createRoadsFromVoronoiEdges(voronoiEdges, roads, m_program);

    // Place buildings within Voronoi cells
    placeBuildingsInVoronoiCells(voronoiCells, gen, m_program, 10.0f);
}

void CityGen::deGenerate()
{
    roads.clear();
    buildings.clear();
}
void CityGen::reGenerate()
{
    deGenerate();
    generate(m_program);
}

void CityGen::render(const glm::mat4& proj, const glm::mat4& view, float m_timer)
{
    for (auto& road : roads)
    {
        road.draw(proj, view, m_timer);
    }
    for (auto& building : buildings)
    {
        building.draw(proj, view, m_timer);
    }
}

*/

/*
#include "CityGen.h"
#include "Building.h"
#include "Road.h"
#include "Cube.h"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
using namespace std;
using namespace glm;

vector<Building> buildings;
vector<Road> roads;

const int gridSize = 100;
const float spacing = 30.0f;
const float citySpacing = 1.1f;
const int sectorSize = 25;
const int skyScraperSector = 5;

const int numSites = 15;    

struct Site {
    float x, z; // Using float for more precise positioning
};

vector<Site> generateVoronoiSites(int numSites, int gridSize, mt19937& gen)
{
    uniform_real_distribution<float> siteDist(-(gridSize * citySpacing), gridSize * citySpacing);
    vector<Site> sites;
    for (int i = 0; i < numSites; ++i) {
        sites.push_back({ siteDist(gen), siteDist(gen) });
    }
    return sites;
}

struct Cell {
    int siteIndex;
};

// Function to assign each grid point to the nearest Voronoi site
vector<vector<int>> assignGridToSites(int gridSize, const vector<Site>& sites, float spacing) {
    vector<vector<int>> grid(gridSize, vector<int>(gridSize, -1));
    float originX = - (gridSize / 2) * spacing;
    float originZ = - (gridSize / 2) * spacing;
    for (int x = 0; x < gridSize; ++x) {
        for (int z = 0; z < gridSize; ++z) {
            float worldX = originX + x * spacing;
            float worldZ = originZ + z * spacing;
            float minDist = numeric_limits<float>::max();
            int closestSite = -1;
            for (int i = 0; i < sites.size(); ++i) {
                float dx = worldX - sites[i].x;
                float dz = worldZ - sites[i].z;
                float dist = dx * dx + dz * dz;
                if (dist < minDist) {
                    minDist = dist;
                    closestSite = i;
                }
            }
            grid[x][z] = closestSite;
        }
    }
    return grid;
}

// Function to extract Voronoi edges based on grid assignments
vector<pair<vec3, vec3>> extractVoronoiEdges(const vector<vector<int>>& grid, int gridSize, float spacing) {
    vector<pair<vec3, vec3>> edges;
    
    float originX = -(gridSize / 2) * spacing;
    float originZ = -(gridSize / 2) * spacing;

    for (int x = 0; x < gridSize - 1; ++x) {
        for (int z = 0; z < gridSize - 1; ++z) {
            
            int current = grid[x][z];
            int right = grid[x + 1][z];
            int up = grid[x][z + 1];
            
            if (current != right) {
                // Vertical edge between (x, z) and (x+1, z)
                vec3 p1(originX + x * spacing, 0.0f, originZ + z * spacing);
                vec3 p2(originX + (x + 1) * spacing, 0.0f, originZ + z * spacing);
                edges.emplace_back(p1, p2);
            }
            if (current != up) {
                // Horizontal edge between (x, z) and (x, z+1)
                vec3 p1(originX + x * spacing, 0.0f, originZ + z * spacing);
                vec3 p2(originX + x * spacing, 0.0f, originZ + (z + 1) * spacing);
                edges.emplace_back(p1, p2);
            }
        }
    }
    return edges;
}

// Function to create roads from Voronoi edges
void createRoadsFromVoronoiEdges(const vector<pair<vec3, vec3>>& edges, vector<Road>& roads, GLuint program) {
    for (const auto& edge : edges) {
        // Depending on your Road implementation, you might need to create multiple Road segments
        // along the edge or represent the entire edge as a single Road.
        // Here, we'll create a Road at the midpoint for simplicity.

        glm::vec3 direction = glm::normalize(edge.first - edge.second);

        vec3 midpoint = (edge.first + edge.second) / 2.0f;
        Road road(midpoint, 0);
        road.init(program);
        roads.push_back(road);

        // Road road1(edge.first,0);
        // road1.init(program);
        // roads.push_back(road1);

        // Road road2(edge.second,0);
        // road2.init(program);
        // roads.push_back(road2);
    }
}

// Function to place buildings around Voronoi sites within their respective cells
void placeBuildingsAroundSites(const vector<Site>& sites, float spacing, mt19937& gen, GLuint program, float scaleHeight) {
    uniform_real_distribution<float> offsetDist(-spacing / 2, spacing / 2);
    uniform_real_distribution<float> heightDist(1.0f, scaleHeight);
    
    for (const auto& site : sites) {
        // Define number of buildings per site
        int numBuildings = 10; // Adjust as needed based on site density
        
        for (int i = 0; i < numBuildings; ++i) {
            float offsetX = offsetDist(gen);
            float offsetZ = offsetDist(gen);
            vec3 buildingPosition = vec3(site.x + offsetX, 0.0f, site.z + offsetZ);
            
            // Determine building height
            float height = heightDist(gen);
            
            Building b(buildingPosition, height);
            b.init(program);
            buildings.push_back(b);
        }
    }
}

void CityGen::generate(GLuint program)
{
    m_program = program;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> scaleDist(0.4f, 2.0f);
    uniform_real_distribution<> scaleHeight(1.0f, 2.0f);
    uniform_real_distribution<> offsetDist(-1.0f, 1.0f);

    vector<Site> sites = generateVoronoiSites(numSites, gridSize, gen);
    
    if (sites.size() >= 4) {
        sites.push_back(sites[0]);
        sites.push_back(sites[1]);
        sites.push_back(sites[2]);
    }
    vector<vector<int>> grid = assignGridToSites(gridSize, sites, spacing);

    vector<pair<vec3, vec3>> voronoiEdges = extractVoronoiEdges(grid, gridSize, spacing);

    createRoadsFromVoronoiEdges(voronoiEdges, roads, m_program);
    
    placeBuildingsAroundSites(sites, spacing, gen, m_program, 10.0f);
    
}

void CityGen::reGenerate()
{
    roads.clear();
    buildings.clear();
    generate(m_program);
}
void CityGen::deGenerate()
{
    roads.clear();
    buildings.clear();
}

void CityGen::render(const glm::mat4& proj, const glm::mat4& view, float m_timer)
{
    // glUseProgram(m_program);

    // glUniformMatrix4fv(glGetUniformLocation(m_program, "projection"), 1, GL_FALSE, glm::value_ptr(proj));
    // glUniformMatrix4fv(glGetUniformLocation(m_program, "view"), 1, GL_FALSE, glm::value_ptr(view));
    // glUniform1f(glGetUniformLocation(m_program, "u_time"), m_timer);

    // glBindVertexArray(buildings[0].getVAO());

    // glDrawArraysInstanced(GL_TRIANGLES, 0, buildings[0].getVertexCount(), instanceTransforms.size());

    for (auto& road : roads)
    {
        road.draw(proj, view, m_timer);
    }
    for (auto& building : buildings)
    {
        building.draw(proj, view, m_timer);
    }
}
*/