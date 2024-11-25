
/*
#include "Building.h"
#include <glm/gtc/type_ptr.hpp>

Building::Building()
    : m_vao(0), m_vbo(0), m_instanceVBO(0), m_program(0) {}

void Building::addInstance(const glm::vec3& position, float scale) {
    m_instanceData.emplace_back(position, scale);
}

void Building::init(GLuint program) {
    m_program = program;
    generateMesh();

    // Generate VAO
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    // Generate VBO for vertex data
    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size() * sizeof(Vertex), m_vertices.data(), GL_STATIC_DRAW);

    // Vertex attributes
    // Position
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, x));
    glEnableVertexAttribArray(0);
    // Color
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (void*)offsetof(Vertex, r));
    glEnableVertexAttribArray(1);

    // Generate VBO for instance data
    glGenBuffers(1, &m_instanceVBO);
    glBindBuffer(GL_ARRAY_BUFFER, m_instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, m_instanceData.size() * sizeof(glm::vec4), m_instanceData.data(), GL_STATIC_DRAW);

    // Instance attributes
    // Position (vec3) and Scale (float) packed into vec4
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1); // Update per instance

    glBindVertexArray(0);
}

void Building::generateMesh() {
    // Define vertices for a cube or building base
    // For simplicity, let's define a unit cube
    GLubyte r = 150, g = 150, b = 150, a = 255;

    // Vertices of a unit cube centered at the origin
    m_vertices = {
        // Front face
        { -0.5f, 0.0f,  0.5f, r, g, b, a },
        {  0.5f, 0.0f,  0.5f, r, g, b, a },
        {  0.5f, 1.0f,  0.5f, r, g, b, a },
        { -0.5f, 1.0f,  0.5f, r, g, b, a },
        // Back face
        { -0.5f, 0.0f, -0.5f, r, g, b, a },
        {  0.5f, 0.0f, -0.5f, r, g, b, a },
        {  0.5f, 1.0f, -0.5f, r, g, b, a },
        { -0.5f, 1.0f, -0.5f, r, g, b, a },
        // ... (other faces)
    };
}

void Building::draw(const glm::mat4& proj, const glm::mat4& view, float m_timer) {
    glUseProgram(m_program);

    // Set uniform variables
    GLint projLoc = glGetUniformLocation(m_program, "projection");
    GLint viewLoc = glGetUniformLocation(m_program, "view");

    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));

    // Update instance VBO if necessary
    glBindBuffer(GL_ARRAY_BUFFER, m_instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, m_instanceData.size() * sizeof(glm::vec4), m_instanceData.data(), GL_STATIC_DRAW);

    // Draw instances
    glBindVertexArray(m_vao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(m_vertices.size()), static_cast<GLsizei>(m_instanceData.size()));
    glBindVertexArray(0);

    glUseProgram(0);
}
*/
/*
// Building.cpp
#include "Building.h"
#include <random>

Building::Building(const std::vector<glm::vec3>& baseVertices)
    : m_baseVertices(baseVertices), m_vao(0), m_vbo(0), m_ebo(0), m_program(0) {
    // Random height between 10 and 50 units
    std::mt19937 gen(static_cast<unsigned int>(time(nullptr)));
    std::uniform_real_distribution<float> dist(10.0f, 50.0f);
    // m_height = dist(gen);
}

void Building::init(GLuint program) {
    m_program = program;
    generateMesh();

    // Create and bind VAO and VBO
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    // Vertex buffer
    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, m_vertices.size() * sizeof(Vertex), m_vertices.data(), GL_STATIC_DRAW);

    // Element buffer
    glGenBuffers(1, &m_ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned int), m_indices.data(), GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(0);

    // Color attribute
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (void*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);

    // Unbind VAO
    glBindVertexArray(0);
}

void Building::generateMesh() {
    // Generate vertices for base and top
    GLubyte r = 150, g = 255, b = 150, a = 255; // Gray color
    size_t n = m_baseVertices.size();

    // Base vertices
    for (const auto& v : m_baseVertices) {
        m_vertices.push_back({ v.x, v.y, v.z, r, g, b, a });
    }

    // Top vertices
    for (const auto& v : m_baseVertices) {
        m_vertices.push_back({ v.x, v.y + m_height, v.z, r, g, b, a });
    }

    // Generate side faces
    for (size_t i = 0; i < n; ++i) {
        unsigned int idx0 = i;
        unsigned int idx1 = (i + 1) % n;
        unsigned int idx2 = i + n;
        unsigned int idx3 = ((i + 1) % n) + n;

        // First triangle
        m_indices.push_back(idx0);
        m_indices.push_back(idx2);
        m_indices.push_back(idx1);

        // Second triangle
        m_indices.push_back(idx1);
        m_indices.push_back(idx2);
        m_indices.push_back(idx3);
    }

    // Optionally, generate top face
    // Triangulate the top face
    for (size_t i = 1; i < n - 1; ++i) {
        m_indices.push_back(n);
        m_indices.push_back(n + i);
        m_indices.push_back(n + i + 1);
    }
}

mat4 Building::getTransformationMatrix() const {
    mat4 model = translate(mat4(1.0f), m_position);
    model = rotate(model, radians(-90.0f), vec3(1.0f, 0.0f, 0.0f));

    // model = rotate(model, radians(m_rotation.x), m_rotationAxis);
    // model = rotate(model, radians(m_rotation.y), m_rotationAxis);
    // model = rotate(model, radians(m_rotation.z), m_rotationAxis);

    model = scale(model, m_scale);
    return model;
}

void Building::draw(const glm::mat4& proj, const glm::mat4& view, float m_timer) {
    glUseProgram(m_program);

    // Set uniform variables
    GLint projLoc = glGetUniformLocation(m_program, "projection");
    GLint viewLoc = glGetUniformLocation(m_program, "view");
    GLint modelLoc = glGetUniformLocation(m_program, "model");

    // glm::mat4 model = glm::mat4(1.0f);
    
    mat4 model = getTransformationMatrix();

    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    // Bind VAO and draw elements
    glBindVertexArray(m_vao);
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_indices.size()), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glUseProgram(0);
}
*/



// /*
#include "Building.h"
#include "Cube.h"
#include "stb_image.h"
#include <vector>
#include <random>
using namespace std;
using namespace glm;

Building::Building() : m_position(vec3(0.0f)), m_scale(vec3(1.0f)), cube(vec3(0.0f), vec3(1.0f), vec3(0.0f), vec4(1.0f)) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> scaleDist(0.4f, 2.0f);
    uniform_real_distribution<> heightDist(1.0f, 3.0f);
    uniform_real_distribution<> offsetDist(-1.0f, 1.0f);

    // Generate random values for scale and position
    float randomHeight = heightDist(gen);
    float randomWidth = scaleDist(gen);
    float randomDepth = scaleDist(gen);

    m_scale = vec3(randomWidth, randomHeight * 3.0f, randomDepth);     // Random scale (width, height, depth)

    // Initialize the cube with these random parameters
    cube = Cube(vec3(0.0f), m_scale, vec3(0.0f), vec4(1.0f));
    cube1 = Cube(vec3(0.0f), m_scale, vec3(0.0f), vec4(1.0f));
}

Building::Building(vec3 pos, float height) : m_position(pos), m_scale(vec3(1.0f)), cube(pos, vec3(1.0f), vec3(0.0f), vec4(1.0f))
{
    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<> heightDist(1.0f, 3.0f);
    float randomHeight = heightDist(gen) * height;//* 2.0f;

    uniform_real_distribution<> scaleDist(1.0f, 2.0f);
    float randomWidth = scaleDist(gen);
    float randomDepth = scaleDist(gen);

    m_scale = vec3(randomWidth, randomDepth, randomHeight);
    m_position = vec3(pos.x, randomHeight / 2.0f , pos.z);
    cube = Cube(m_position, m_scale, vec3(0.0f), vec4(1.0f));
    
    uniform_real_distribution<> secondOffsett(-1.0f, 1.0f);
    cube1 = Cube(m_position + vec3(secondOffsett(gen),randomHeight * 0.5f, secondOffsett(gen)), m_scale * 0.5f, vec3(0.0f), vec4(1.0f));
    
}

void Building::setPosition(const glm::vec3& pos)
{
    m_position = pos;
    cube.setPosition(m_position);
    cube1.setPosition(m_position);
}

void Building::init(GLuint program)
{
    cube.init(program);
    cube1.init(program);
}
void Building::draw(const glm::mat4& proj, const glm::mat4& view, float m_timer)
{
    cube.draw(proj,view,m_timer);
}

void randomBuilding()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> scaleDist(0.4f, 2.0f);
    uniform_real_distribution<> scaleHeight(1.0f, 2.0f);
    uniform_real_distribution<> offsetDist(-1.0f, 1.0f);

    uniform_real_distribution<> offsetRoad(3.0f, 8.0f);

    // int xRoadOffset = offsetRoad(gen);
    // int zRoadOffset = offsetRoad(gen);

    // float height = scaleHeight(gen) * 3.0f;

    vec3 scale(scaleDist(gen) * 1.25f,
    scaleDist(gen) * 1.25f,
    0.01f);//height);

    // x = Cube(vec3(0.0f), scale, vec3(0.0f, 0.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));

    // cube.setPosition(position);

    // cube.init(m_program);

    // buildings.push_back(cube);    
}
// */