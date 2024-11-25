/* 3rd
#pragma once
#include <vector>
#include "Vertex.h"
#include <glm/glm.hpp>
#include <GL/glew.h>
using namespace std;
using namespace glm;

class Building {
public:
    Building();
    void init(GLuint program);
    void draw(const glm::mat4& proj, const glm::mat4& view, float m_timer);

    void addInstance(const glm::vec3& position, float scale);

private:
    vector<Vertex> m_vertices;
    vector<glm::vec4> m_instanceData; // Holds position and scale
    GLuint m_vao, m_vbo, m_instanceVBO;
    GLuint m_program;

    void generateMesh();
};
*/

/* 2nd
#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Cube.h"
using namespace std;
using namespace glm;

class Building {
private:
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_ebo;
    GLuint m_program;
    float m_height = 5.0f;
    vector<glm::vec3> m_baseVertices;
    vector<Vertex> m_vertices;
    vector<int> m_indices;

    Cube cube;
    Cube cube1;
    glm::vec3 m_position;
    glm::vec3 m_scale;
public:

    Building();
    Building(const vector<vec3>& baseVertices);

    void init(GLuint program);
    void draw(const glm::mat4& proj, const glm::mat4& view, float m_timer);

    void generateMesh();
    mat4 getTransformationMatrix() const;
    
    void randomBuilding();
    void setPosition(const glm::vec3& pos);
};
*/

// /* 1st
#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Cube.h"
using namespace std;

class Building {
private:
    Cube cube;
    Cube cube1;
    glm::vec3 m_position;
    glm::vec3 m_scale;
public:

    Building();
    Building(glm::vec3 pos, float height);

    void init(GLuint program);
    void draw(const glm::mat4& proj, const glm::mat4& view, float m_timer);

    void randomBuilding();
    void setPosition(const glm::vec3& pos);
};
// */