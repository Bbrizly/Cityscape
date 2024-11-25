#include "CityGen.h"
using namespace std;
using namespace glm;
/*
Cube::Cube()
    : m_position(vec3(0.0f)), m_scale(vec3(1.0f)), m_rotation(vec3(0.0f)), m_color(vec4(1.0f))
{
    createQuads();
}

// Parameterized constructor
Cube::Cube(const vec3& position, const vec3& scale, const vec3& rotation, const vec4& color)
    : m_position(position), m_scale(scale), m_rotation(rotation), m_color(color)
{
    createQuads();
}

// Initialize all quads
void Cube::init(GLuint program) {
    for (auto& quad : m_quads) {
        quad.init(program);
    }
}

// Render the cube
void Cube::draw(const mat4& proj, const mat4& view, float m_timer, GLuint program) {
    // Calculate the model matrix
    mat4 model = getModelMatrix();

    // Use the shader program
    glUseProgram(program);

    // Set uniform matrices
    GLint projLoc = glGetUniformLocation(program, "projection");
    GLint viewLoc = glGetUniformLocation(program, "view");
    GLint worldLoc = glGetUniformLocation(program, "world");
    GLint timeLoc  = glGetUniformLocation(program, "u_time");

    glUniformMatrix4fv(projLoc, 1, GL_FALSE, value_ptr(proj));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, value_ptr(view));
    glUniformMatrix4fv(worldLoc, 1, GL_FALSE, value_ptr(model));
    glUniform1f(timeLoc, m_timer);

    // Draw each quad
    for (auto& quad : m_quads) {
        quad.draw(mat4(1.0f), program); // Pass identity matrix since transformations are handled in the model matrix
    }
}

// Setters for transformation properties
void Cube::setPosition(const vec3& pos) { m_position = pos; }
void Cube::setScale(const vec3& scl) { m_scale = scl; }
void Cube::setRotation(const vec3& rot) { m_rotation = rot; }

// Set the color of the cube and update all quads
void Cube::setColor(const vec4& color) {
    m_color = color;
    for (auto& quad : m_quads) {
        quad.setColor(color);
    }
}

// Calculate the model matrix based on position, scale, and rotation
mat4 Cube::getModelMatrix() const {
    mat4 model = mat4(1.0f);

    // Apply translations
    model = translate(model, m_position);

    // Apply rotations (Euler angles in degrees)
    model = rotate(model, radians(m_rotation.x), vec3(1.0f, 0.0f, 0.0f));
    model = rotate(model, radians(m_rotation.y), vec3(0.0f, 1.0f, 0.0f));
    model = rotate(model, radians(m_rotation.z), vec3(0.0f, 0.0f, 1.0f));

    // Apply scaling
    model = scale(model, m_scale);

    return model;
}

// Helper function to create six quads for the cube
void Cube::createQuads() {
    // Clear existing quads if any
    m_quads.clear();

    // Six faces of the cube
    // Front (+Z)
    m_quads.emplace_back(Quad(m_color));
    // Back (-Z)
    m_quads.emplace_back(Quad(m_color));
    // Left (-X)
    m_quads.emplace_back(Quad(m_color));
    // Right (+X)
    m_quads.emplace_back(Quad(m_color));
    // Top (+Y)
    m_quads.emplace_back(Quad(m_color));
    // Bottom (-Y)
    m_quads.emplace_back(Quad(m_color));
}
*/
vector<Vertex> m_vertices = {
    { -0.5f, -0.5f,  0.5f, 255, 0, 0, 255 },
    { -0.5f,  0.5f,  0.5f, 255, 0, 0, 255 },
    {  0.5f,  0.5f,  0.5f, 255, 0, 0, 255 },
    {  0.5f,  0.5f,  0.5f, 255, 0, 0, 255 },
    {  0.5f, -0.5f,  0.5f, 255, 0, 0, 255 },
    { -0.5f, -0.5f,  0.5f, 255, 0, 0, 255 },

    {  0.5f,  0.5f, -0.5f, 128, 0, 0, 255 },
    { -0.5f,  0.5f, -0.5f, 128, 0, 0, 255 },
    { -0.5f, -0.5f, -0.5f, 128, 0, 0, 255 },
    { -0.5f, -0.5f, -0.5f, 128, 0, 0, 255 },
    {  0.5f, -0.5f, -0.5f, 128, 0, 0, 255 },
    {  0.5f,  0.5f, -0.5f, 128, 0, 0, 255 },

    { -0.5f,  0.5f, -0.5f, 0, 255, 0, 255 },
    { -0.5f,  0.5f,  0.5f, 0, 255, 0, 255 },
    { -0.5f, -0.5f,  0.5f, 0, 255, 0, 255 },
    { -0.5f, -0.5f,  0.5f, 0, 255, 0, 255 },
    { -0.5f, -0.5f, -0.5f, 0, 255, 0, 255 },
    { -0.5f,  0.5f, -0.5f, 0, 255, 0, 255 },

    {  0.5f,  0.5f,  0.5f, 0, 128, 0, 255 },
    {  0.5f,  0.5f, -0.5f, 0, 128, 0, 255 },
    {  0.5f, -0.5f, -0.5f, 0, 128, 0, 255 },
    {  0.5f, -0.5f, -0.5f, 0, 128, 0, 255 },
    {  0.5f, -0.5f,  0.5f, 0, 128, 0, 255 },
    {  0.5f,  0.5f,  0.5f, 0, 128, 0, 255 },

    { -0.5f,  0.5f,  0.5f, 0, 0, 255, 255 },
    { -0.5f,  0.5f, -0.5f, 0, 0, 255, 255 },
    {  0.5f,  0.5f, -0.5f, 0, 0, 255, 255 },
    {  0.5f,  0.5f, -0.5f, 0, 0, 255, 255 },
    {  0.5f,  0.5f,  0.5f, 0, 0, 255, 255 },
    { -0.5f,  0.5f,  0.5f, 0, 0, 255, 255 },

    { -0.5f, -0.5f,  0.5f, 0, 0, 128, 255 },
    {  0.5f, -0.5f,  0.5f, 0, 0, 128, 255 },
    {  0.5f, -0.5f, -0.5f, 0, 0, 128, 255 },
    {  0.5f, -0.5f, -0.5f, 0, 0, 128, 255 },
    { -0.5f, -0.5f, -0.5f, 0, 0, 128, 255 },
    { -0.5f, -0.5f,  0.5f, 0, 0, 128, 255 },
};

Cube::Cube()
    : m_position(vec3(0.0f)), m_scale(vec3(1.0f)), m_rotation(vec3(1.0f)), m_color(vec4(1.0f)) {}

Cube::Cube(const vec3& position, const vec3& scale, const vec3& rotation, const vec4& color)
    : m_position(position), m_scale(scale), m_rotation(rotation), m_color(color) {}

// Calculates and returns the transformation matrix for the cube, combining translation, rotation, and scaling
mat4 Cube::getTransformationMatrix() const {
    mat4 model = translate(mat4(1.0f), m_position);
    model = rotate(model, radians(-90.0f), vec3(1.0f, 0.0f, 0.0f));

    // model = rotate(model, radians(m_rotation.x), m_rotationAxis);
    // model = rotate(model, radians(m_rotation.y), m_rotationAxis);
    // model = rotate(model, radians(m_rotation.z), m_rotationAxis);

    model = scale(model, m_scale);
    return model;
}

// Rotates the cube based on rotation axis and speed, adjusted by deltaTime
// void Cube::rotating(float deltaTime) 
// {
//     m_rotation += m_rotationAxis * m_rotationSpeed * deltaTime;
// }

void Cube::init(GLuint program) {
    m_program = program;
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * m_vertices.size(), m_vertices.data(), GL_DYNAMIC_DRAW);

    int posAttr = glGetAttribLocation(program, "a_position");
    int colorAttr = glGetAttribLocation(program, "a_color");

    glVertexAttribPointer(posAttr, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(posAttr);

    glVertexAttribPointer(colorAttr, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (void*)(sizeof(GLfloat) * 3));
    glEnableVertexAttribArray(colorAttr);
}

void Cube::draw(const mat4& proj, const mat4& view, float m_timer) {
    glUseProgram(m_program);
    glBindVertexArray(m_vao);

    mat4 model = getTransformationMatrix();
    // GLint program;
    // glGetIntegerv(GL_CURRENT_PROGRAM,&program);

    glUniformMatrix4fv(glGetUniformLocation(m_program, "projection"), 1, GL_FALSE, value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "view"), 1, GL_FALSE, value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "world"), 1, GL_FALSE, value_ptr(model));
    glUniform1f(glGetUniformLocation(m_program, "u_time"), m_timer);

    glDrawArrays(GL_TRIANGLES, 0, m_vertices.size());

    glBindVertexArray(0);
    glUseProgram(0);
}

// Setters for cube transformation properties: position, rotation axis, rotation and scale
void Cube::setPosition(const vec3& pos) { m_position = pos; }
// void Cube::setRotationAxis(const vec3& rot) { m_rotationAxis = rot; }
void Cube::setRotation(const vec3& rot) { m_rotation = rot; }
void Cube::setScale(const vec3& scl) { m_scale = scl; }

// // cube color getter
// vec4 Cube::getColor() { return m_color; }
// vec3 Cube::getScale() { return m_scale; }

// cube color setter
void Cube::setColor(const vec4& col) 
{
    m_color = col;
    
    for (auto& vertex : m_vertices) {
        vertex.r = static_cast<GLubyte>(col.r * 255);
        vertex.g = static_cast<GLubyte>(col.g * 255);
        vertex.b = static_cast<GLubyte>(col.b * 255);
        vertex.a = static_cast<GLubyte>(col.a * 255);
    }

    // Re-upload the vertex buffer data to the GPU
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertex) * m_vertices.size(), m_vertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// Returns pointer to the cube vertices array
vector<Vertex> Cube::getVertices() { return m_vertices; }

// Returns the total number  of vertices in the cube
int Cube::getVertexCount() { return sizeof(m_vertices) / sizeof(Vertex); }