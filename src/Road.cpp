#include "Road.h"
#include "Cube.h"
#include <vector>
#include <random>
using namespace std;
using namespace glm;

Road::Road() : m_position(vec3(0.0f)), m_scale(vec3(1.0f)), cube(vec3(0.0f), vec3(1.0f), vec3(0.0f), vec4(1.0f)) {
    m_scale = vec3(2.0f,2.0f,0.0f);
    m_position = vec3(0.0f);
    cube = Cube(m_position, m_scale, vec3(0.0f), vec4(1.0f));
}

Road::Road(vec3 pos, int interRightLeft) : intersectionRightLeft(interRightLeft), m_position(pos), m_scale(vec3(1.0f)), cube(pos, vec3(1.0f), vec3(0.0f), vec4(1.0f))
{
    //interRightLeft = 
    //0 for intersection,
    //1 for Right
    //2 for Left
    float roadSize = 3.0f;
    m_scale = vec3(roadSize, roadSize, 0.0f);
    m_position = vec3(pos.x, pos.y, pos.z);
    cube = Cube(m_position, m_scale, vec3(0.0f), vec4(1.0f));
}

void Road::test()
{
    // return;

    if(intersectionRightLeft == 0)          //intersection
    {
        cube.setColor(vec4(255,255,0,255));
    }
    else if(intersectionRightLeft == 1)     // Along x
    {
        cube.setColor(vec4(51,255,255,255));
    }
    else if(intersectionRightLeft == 2)// Along y
    {
        cube.setColor(vec4(128,128,128,255));
    }
}

void Road::setPosition(const glm::vec3& pos)
{
    m_position = pos;
    cube.setPosition(m_position);
}

void Road::init(GLuint program)
{
    // cube.setColor(vec4(128,128,128,255));
    cube.init(program);
}
void Road::draw(const glm::mat4& proj, const glm::mat4& view, float m_timer)
{   
    //
    GLint program;
    glGetIntegerv(GL_CURRENT_PROGRAM,&program);

    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, m_position);
    model = glm::rotate(model, glm::radians(rotationAngle), glm::vec3(0.0f, 1.0f, 0.0f));
    model = glm::scale(model, m_scale);

    // Assuming your shader has a 'model' matrix uniform
    glUniformMatrix4fv(glGetUniformLocation(program, "model"), 1, GL_FALSE, glm::value_ptr(model));
    //
    
    cube.draw(proj,view,m_timer);
}