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