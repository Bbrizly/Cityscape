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