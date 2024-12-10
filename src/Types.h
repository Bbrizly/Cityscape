#pragma once
#include <GL/glew.h>
#include <vector>
using namespace std;

enum District{
    Industrial,
    Commercial,
    Residential
};

struct Point { double x, y; };
struct Line { double a, b, c; };

struct Building {
    vector<Point> polygons;
    int textureLayer;
    float height;
    int district;
    GLubyte r; GLubyte g; GLubyte b;
    bool isSpecial;
};
