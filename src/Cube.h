#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "../wolf/wolf.h"
#include "Vertex.h"
#include "Quad.h"
using namespace std;

class Cube {
private:
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_program;

    glm::vec3 m_position;
    glm::vec3 m_scale;
    glm::vec3 m_rotation;
    glm::vec4 m_color;

    void createQuads();
    std::vector<Quad> m_quads;

    glm::mat4 getModelMatrix() const;

public:
    Cube();
    Cube(const glm::vec3& position, const glm::vec3& scale, const glm::vec3& rotation, const glm::vec4& color);

    glm::mat4 getTransformationMatrix() const;

    void init(GLuint program);

    void draw(const glm::mat4& proj, const glm::mat4& view, float m_timer);

    void setPosition(const glm::vec3& pos);
    void setScale(const glm::vec3& scl);
    void setRotation(const glm::vec3& rot);
    void setColor(const glm::vec4& col);

    vector<Vertex> getVertices();
    int getVertexCount();
};
