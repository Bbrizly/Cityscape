#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Vertex.h"
using namespace glm;

class Quad {
private:
    std::vector<Vertex> m_vertices;
    GLuint m_vao, m_vbo;
    glm::vec4 m_color;
public:
    Quad(const glm::vec4& color = glm::vec4(1.0f));
    void init(GLuint program);
    void draw(const glm::mat4& model, GLuint program);
    
    void setColor(const glm::vec4& color);
};
