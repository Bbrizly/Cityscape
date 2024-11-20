#include "Quad.h"
#include <glm/gtc/type_ptr.hpp>

// Constructor implementation
Quad::Quad(const glm::vec4& color)
    : m_color(color)
{
    // Define a standard quad centered at the origin in the X-Y plane
    // Z-coordinate is 0 for all vertices
    // Two triangles: (v0, v1, v2) and (v0, v2, v3)

    Vertex v0 = { -0.5f, -0.5f, 0.0f,
                 static_cast<GLubyte>(color.r * 255),
                 static_cast<GLubyte>(color.g * 255),
                 static_cast<GLubyte>(color.b * 255),
                 static_cast<GLubyte>(color.a * 255) };

    Vertex v1 = { -0.5f,  0.5f, 0.0f,
                 static_cast<GLubyte>(color.r * 255),
                 static_cast<GLubyte>(color.g * 255),
                 static_cast<GLubyte>(color.b * 255),
                 static_cast<GLubyte>(color.a * 255) };

    Vertex v2 = {  0.5f,  0.5f, 0.0f,
                 static_cast<GLubyte>(color.r * 255),
                 static_cast<GLubyte>(color.g * 255),
                 static_cast<GLubyte>(color.b * 255),
                 static_cast<GLubyte>(color.a * 255) };

    Vertex v3 = {  0.5f, -0.5f, 0.0f,
                 static_cast<GLubyte>(color.r * 255),
                 static_cast<GLubyte>(color.g * 255),
                 static_cast<GLubyte>(color.b * 255),
                 static_cast<GLubyte>(color.a * 255) };

    // First triangle
    m_vertices.push_back(v0);
    m_vertices.push_back(v1);
    m_vertices.push_back(v2);

    // Second triangle
    m_vertices.push_back(v0);
    m_vertices.push_back(v2);
    m_vertices.push_back(v3);
}

// Initialize OpenGL buffers
void Quad::init(GLuint program) {
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * m_vertices.size(), m_vertices.data(), GL_STATIC_DRAW);

    // Define attribute locations
    // Assuming the shader uses location 0 for position and 1 for color
    GLint posAttr = 0;    // Position attribute location
    GLint colorAttr = 1;  // Color attribute location

    // Position attribute
    glVertexAttribPointer(posAttr, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(posAttr);

    // Color attribute
    glVertexAttribPointer(colorAttr, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (void*)(offsetof(Vertex, r)));
    glEnableVertexAttribArray(colorAttr);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

// Render the quad
void Quad::draw(const glm::mat4& model, GLuint program) {
    glUseProgram(program);
    
    // Set the model matrix uniform
    GLint modelLoc = glGetUniformLocation(program, "world");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    glBindVertexArray(m_vao);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(m_vertices.size()));
    glBindVertexArray(0);
}

// Optional: Update color dynamically
void Quad::setColor(const glm::vec4& color) {
    for (auto& vertex : m_vertices) {
        vertex.r = static_cast<GLubyte>(color.r * 255);
        vertex.g = static_cast<GLubyte>(color.g * 255);
        vertex.b = static_cast<GLubyte>(color.b * 255);
        vertex.a = static_cast<GLubyte>(color.a * 255);
    }
    
    // Update the vertex buffer
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertex) * m_vertices.size(), m_vertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}