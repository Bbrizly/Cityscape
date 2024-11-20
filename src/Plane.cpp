#include "Plane.h"
using namespace std;

vector<Vertex> generatePlaneVertices(int subdivisions) {
    vector<Vertex> vertices;
    float step = 1.0f / subdivisions;

    for (int y = 0; y < subdivisions; ++y) {
        for (int x = 0; x < subdivisions; ++x) {
            
            float x0 = x * step - 0.5f;
            float y0 = y * step - 0.5f;
            float x1 = (x + 1) * step - 0.5f;
            float y1 = (y + 1) * step - 0.5f;

            vertices.push_back({x0, y0, 0.0f, 0, 0, 0, 255});
            vertices.push_back({x1, y0, 0.0f, 0, 0, 0, 255});
            vertices.push_back({x0, y1, 0.0f, 0, 0, 0, 255});

            vertices.push_back({x1, y0, 0.0f, 0, 0, 0, 255});
            vertices.push_back({x1, y1, 0.0f, 0, 0, 0, 255});
            vertices.push_back({x0, y1, 0.0f, 0, 0, 0, 255});
        }
    }
    return vertices;
}

glm::mat4 Plane::getTransformationMatrix() const {
    glm::mat4 model = glm::translate(glm::mat4(1.0f), m_position);
    model = glm::rotate(model, glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    model = glm::scale(model, m_scale); 
    return model;
}

Plane::Plane()
    : m_position(0.0f), m_scale(1.0f), m_color(1.0f), m_subdivisions(1) {
    m_vertices = generatePlaneVertices(m_subdivisions);}
Plane::Plane(const glm::vec3& position, const glm::vec3& scale, const glm::vec4& color, int subdivisions)
    : m_position(position), m_scale(scale), m_color(color), m_subdivisions(subdivisions) {
    m_vertices = generatePlaneVertices(m_subdivisions);}

void Plane::init(GLuint program) {
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

void Plane::draw(const glm::mat4& proj, const glm::mat4& view, float m_timer) const {
    glBindVertexArray(m_vao);

    glm::mat4 model = getTransformationMatrix();
    GLint program;
    glGetIntegerv(GL_CURRENT_PROGRAM, &program);

    glUniformMatrix4fv(glGetUniformLocation(program, "projection"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(program, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(program, "world"), 1, GL_FALSE, glm::value_ptr(model));
    glUniform1f(glGetUniformLocation(program, "u_time"), m_timer);

    glDrawArrays(GL_TRIANGLES, 0, m_vertices.size());
}

void Plane::setSubdivisions(GLuint program, int subs) {
    m_subdivisions = subs;  //set new subdivisions
    m_vertices = generatePlaneVertices(m_subdivisions); //create triangless

    if (program != 0) {
        init(program); //reinitialize the whole thing
    }
}

void Plane::PressOne()
{
    m_currentSpeedIndex = (m_currentSpeedIndex + 1) % 3;
    
    glUseProgram(m_program);
    glUniform1f(glGetUniformLocation(m_program, "speed"), m_speeds[m_currentSpeedIndex]);  
}
void Plane::incrementFrequency()
{
    m_currentFrequencyIndex = (m_currentFrequencyIndex + 1) % 3;
    
    glUseProgram(m_program);
    glUniform1f(glGetUniformLocation(m_program, "frequency"), m_frequencies[m_currentFrequencyIndex]);  
}
void Plane::incrementAmplitude()
{
    // glGetUniformLocation(m_program, "speed") = ;
    m_currentAmplitudeIndex = (m_currentAmplitudeIndex + 1) % 3;

    glUseProgram(m_program);
    glUniform1f(glGetUniformLocation(m_program, "amplitude"), m_amplitudes[m_currentAmplitudeIndex]);    
}
void Plane::incrementColor()
{
    currentColor = !currentColor ? true : false;
    
    glUseProgram(m_program);
    glUniform1f(glGetUniformLocation(m_program, "chooseColor"), currentColor);
}
void Plane::incrementSubdivisions()
{
    m_currentSubdivisionsIndex = (m_currentSubdivisionsIndex + 1) % 3;

    setSubdivisions(m_program, m_subs[m_currentSubdivisionsIndex]);
}

void Plane::setPosition(const glm::vec3& pos) { m_position = pos; }
void Plane::setScale(const glm::vec3& scl) { m_scale = scl; }
void Plane::setColor(const glm::vec4& col) { m_color = col; }

glm::vec4 Plane::getColor() const { return m_color; }
glm::vec3 Plane::getScale() const { return m_scale; }

const Vertex* Plane::getVertices() const { return m_vertices.data(); }

int Plane::getVertexCount() const { return m_vertices.size(); }
