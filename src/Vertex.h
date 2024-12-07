#pragma once
#include <GL/glew.h>

#pragma pack(push, 1)
struct Vertex {
    GLfloat x, y, z;
    GLubyte r, g, b, a;
    GLfloat u, v;
    GLfloat nx, ny, nz;
};
#pragma pack(pop)