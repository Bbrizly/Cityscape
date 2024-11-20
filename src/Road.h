#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "Cube.h"
using namespace std;

class Road {
private:
    Cube cube;
    glm::vec3 m_position;
    glm::vec3 m_scale;
    int intersectionRightLeft = 0;
    float rotationAngle;
    
public:

    Road();
    Road(glm::vec3 pos, int interRightLeft);

    void init(GLuint program);
    void draw(const glm::mat4& proj, const glm::mat4& view, float m_timer);
    void test();
    void setPosition(const glm::vec3& pos);
};
