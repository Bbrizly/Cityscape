/*
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <fstream>
#include <sstream>

// Define structures
struct Point {
    double x, y;
};

struct Line {
    double a, b, c; // Line equation: ax + by + c = 0
};

// Parameters
const int numSites = 15;
const double minX = -500.0;
const double maxX = 500.0;
const double minY = -500.0;
const double maxY = 500.0;
const unsigned int seed = 40; // Seed for random number generator
const double epsilon = 1e-9;  // Small value for numerical precision

// Function to generate random sites
std::vector<Point> generateSites(int numSites, unsigned int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> distX(minX + 10, maxX - 10);
    std::uniform_real_distribution<double> distY(minY + 10, maxY - 10);

    std::vector<Point> sites;
    for (int i = 0; i < numSites; ++i) {
        sites.push_back({ distX(gen), distY(gen) });
    }
    return sites;
}

// Function to compute the perpendicular bisector between two points
Line perpendicularBisector(const Point& p1, const Point& p2) {
    double mx = (p1.x + p2.x) / 2.0;
    double my = (p1.y + p2.y) / 2.0;

    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    double a = -dy;
    double b = dx;
    double c = -(a * mx + b * my);

    return { a, b, c };
}

// Function to evaluate which side of the line a point is on
double evaluate(const Line& l, const Point& p) {
    return l.a * p.x + l.b * p.y + l.c;
}

// Function to compute intersection point between a line segment and a line
bool lineSegmentLineIntersection(const Point& p1, const Point& p2, const Line& l, Point& intersection) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    double denominator = l.a * dx + l.b * dy;
    if (std::abs(denominator) < epsilon) {
        return false; // Parallel
    }

    double t = -(l.a * p1.x + l.b * p1.y + l.c) / denominator;
    if (t < -epsilon || t > 1.0 + epsilon) {
        return false; // Outside segment
    }

    intersection.x = p1.x + t * dx;
    intersection.y = p1.y + t * dy;
    return true;
}

// Function to clip polygon with a half-plane
std::vector<Point> clipPolygon(const std::vector<Point>& polygon, const Line& l, const Point& site) {
    std::vector<Point> newPolygon;
    if (polygon.empty()) return newPolygon;

    double siteSide = evaluate(l, site);
    bool keepPositive = siteSide >= 0;

    Point prev = polygon.back();
    double prevEval = evaluate(l, prev);
    bool prevInside = keepPositive ? (prevEval >= 0) : (prevEval < 0);

    for (const Point& curr : polygon) {
        double currEval = evaluate(l, curr);
        bool currInside = keepPositive ? (currEval >= 0) : (currEval < 0);

        if (currInside) {
            if (!prevInside) {
                Point intersection;
                if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                    newPolygon.push_back(intersection);
                }
            }
            newPolygon.push_back(curr);
        } else if (prevInside) {
            Point intersection;
            if (lineSegmentLineIntersection(prev, curr, l, intersection)) {
                newPolygon.push_back(intersection);
            }
        }

        prev = curr;
        prevInside = currInside;
    }

    return newPolygon;
}

int main() {
    // Generate random sites
    std::vector<Point> sites = generateSites(numSites, seed);

    // Generate SVG file
    std::ofstream svgFile("voronoi.svg");
    svgFile << "<svg xmlns='http://www.w3.org/2000/svg' width='800' height='800' viewBox='"
            << minX << " " << minY << " " << (maxX - minX) << " " << (maxY - minY) << "'>\n";

    svgFile << "<rect x='" << minX << "' y='" << minY << "' width='" << (maxX - minX)
            << "' height='" << (maxY - minY) << "' fill='white' stroke='black'/>\n";

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> colorDist(0.2, 0.8);

    for (int i = 0; i < numSites; ++i) {
        const Point& site = sites[i];
        std::vector<Point> cell = {
            { minX, minY },
            { maxX, minY },
            { maxX, maxY },
            { minX, maxY }
        };

        for (int j = 0; j < numSites; ++j) {
            if (i == j) continue;

            const Point& otherSite = sites[j];
            Line bisector = perpendicularBisector(site, otherSite);

            cell = clipPolygon(cell, bisector, site);
            if (cell.empty()) break;
        }

        if (!cell.empty()) {
            double r = colorDist(gen) * 255;
            double g = colorDist(gen) * 255;
            double b = colorDist(gen) * 255;

            svgFile << "<polygon points='";
            for (const Point& p : cell) {
                svgFile << p.x << "," << p.y << " ";
            }
            svgFile << "' fill='rgb(" << r << "," << g << "," << b << ")' stroke='black'/>\n";
        }
    }

    for (const Point& site : sites) {
        svgFile << "<circle cx='" << site.x << "' cy='" << site.y
                << "' r='5' fill='red'/>\n";
    }

    svgFile << "</svg>\n";
    svgFile.close();
    std::cout << "Voronoi diagram saved to 'voronoi.svg'.\n";
    return 0;
}


*/
#include <stdio.h>
#include <iostream>
#include <glm/glm.hpp>
#include "../wolf/wolf.h"
#include "../samplefw/SampleRunner.h"
#include "Program.h"

class Week2: public wolf::App
{
private:
    SampleRunner m_sampleRunner;
    Program* program;
public:
    Week2() : App("Citscape - Bassam")
    {
        program = new Program(this);
        m_sampleRunner.addSample(new Program(this));

        
    }

    ~Week2()
    {
    }

    void update(float dt) override
    {
        m_sampleRunner.update(dt);
        
        if(isKeyJustDown('1') || isKeyJustDown('r'))
        {
            program->PressOne();
        }
        else if(isKeyJustDown('2') || isKeyJustDown('d'))
        {
            program->PressTwo();
        }
    }

    void render() override
    {
        m_sampleRunner.render(m_width, m_height);
    }
};

int main(int, char**) {
    Week2 week2;
    week2.run();
}