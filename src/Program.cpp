#include "Program.h"
#include "CityGen.h"
#include <vector>
#include <random>
using namespace std;
using namespace glm;

Program::~Program()
{
	glDeleteVertexArrays(1, &m_vao);
	glDeleteBuffers(1, &m_vbo);
	glDeleteProgram(cubeShader);
    delete m_pOrbitCam;
    delete m_pCamera;
}

// Cube floorr(vec3(0),vec3(50.0f,50.0f,50.0f),vec3(0),vec4(0));

CityGen city;
void Program::init()
{
    if(!cubeShader || !otherShader)
    {
        glEnable(GL_DEPTH_TEST);
        // // glEnable(GL_CULL_FACE);
        // glDisable(GL_CULL_FACE);
        glDepthFunc(GL_LESS);

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        
        // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        m_pProgram = wolf::ProgramManager::CreateProgram("data/floor.vsh", "data/floor.fsh");
        
        // m_program = wolf::LoadShaders("data/planeAnimation.vsh", "data/planeAnimation.fsh");
        cubeShader = wolf::LoadShaders("data/cube.vsh", "data/cube.fsh");
        otherShader = wolf::LoadShaders("data/floor.vsh", "data/floor.fsh");
        city.generate(cubeShader);


        // m_pCamera = new FirstPersonCamera(m_pApp);

        m_pOrbitCam = new OrbitCamera(m_pApp);
        m_pOrbitCam->focusOn(vec3(-5.0f, -5.0f, -5.0f), vec3(0.0f, 100.0f, 0.0f));
    }
}

void Program::update(float x) 
{
	m_timer += x;
    
    // if(m_pCamera) {m_pCamera->update(x);}
    if (m_pOrbitCam) {m_pOrbitCam->update(x);}
}

void Program::render(int width, int height)
{
    // glClearColor(0.3f, 0.3f, 0.3f, 1.0);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    // glClearColor(255.0f, 255.0f, 255.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(cubeShader);

    // mat4 mProj = m_pCamera->getProjMatrix(width, height);
    // mat4 mView = m_pCamera->getViewMatrix();
    mat4 mProj = m_pOrbitCam->getProjMatrix(width, height);
    mat4 mView = m_pOrbitCam->getViewMatrix();

    // glEnable(GL_DEBUG_OUTPUT);

    city.render(mProj,mView,m_timer);
    // floorr.draw(mProj,mView,m_timer);
}

void Program::PressOne()
{
    city.reGenerate();
}

void Program::PressTwo()
{
    city.deGenerate();
}