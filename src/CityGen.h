#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Cube.h"
using namespace std;

class CityGen {
private:
    GLuint m_vao;
    GLuint m_vbo;
    GLuint m_program = 0;

public:

    // struct Site {float x, z;};

    // void computeCatmullRomSpline(const std::vector<Site>& sites, int pointsPerSegment);
    // void generateRoadMesh(const std::vector<glm::vec3>& splinePoints, float roadWidth);
    

    void generate(GLuint program);
    void reGenerate();
    void render(const glm::mat4& proj, const glm::mat4& view, float m_timer);
};
