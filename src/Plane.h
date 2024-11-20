#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "../wolf/wolf.h"
#include "Vertex.h"

// struct Vertex {
//     GLfloat x, y, z;
//     GLubyte r, g, b, a;
// };

class Plane {
private:
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_program; // Store the shader program to use later for modifying uniforms

    
    glm::vec3 m_position;
    glm::vec3 m_scale;
    glm::vec4 m_color;
    int m_subdivisions;

    std::vector<Vertex> m_vertices;

public:
    Plane();
    Plane(const glm::vec3& position, const glm::vec3& scale, const glm::vec4& color, int subdivisions);

    void init(GLuint program);
    void draw(const glm::mat4& proj, const glm::mat4& view, float m_timer) const;

    glm::mat4 getTransformationMatrix() const;
    
    int m_currentSpeedIndex = 0;
    void PressOne();
    float m_speeds[3] = { 5.0f, 8.0f, 3.0f };

    int m_currentFrequencyIndex = 0;
    void incrementFrequency();
    float m_frequencies[3] = { 5.0f, 8.0f, 2.0f };

    int m_currentAmplitudeIndex = 0;
    void incrementAmplitude();
    float m_amplitudes[3] = { 0.5f, 0.8f, 0.2f };

    void incrementColor();
    bool currentColor = true;

    int m_currentSubdivisionsIndex = 0;
    void incrementSubdivisions();
    float m_subs[3] = { 2.0f, 5.0f, 8.0f };

    // Setters for plane properties
    void setPosition(const glm::vec3& pos);
    void setScale(const glm::vec3& scl);
    void setColor(const glm::vec4& col);
    void setSubdivisions(GLuint program, int subs);

    // Getters for plane properties
    glm::vec4 getColor() const;
    glm::vec3 getScale() const;

    // Accessor for vertices and vertex count
    const Vertex* getVertices() const;
    int getVertexCount() const;
};
