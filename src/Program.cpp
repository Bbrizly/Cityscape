#include "Program.h"
#include "CityGen.h"
#include <vector>
#include <random>
using namespace std;
using namespace glm;

Program::~Program()
{
    delete m_pOrbitCam;
    delete m_pCamera;
}


CityGen city;
void Program::init()
{
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    cubeShader = wolf::LoadShaders("data/cube.vsh", "data/cube.fsh");
    city.generate(cubeShader);


    m_pCamera = new FirstPersonCamera(m_pApp);

    // m_pOrbitCam = new OrbitCamera(m_pApp);
    // m_pOrbitCam->focusOn(vec3(-5.0f, -5.0f, -5.0f), vec3(0.0f, 100.0f, 0.0f));
    
}

void Program::update(float x) 
{
	m_timer += x;
    
    if(m_pCamera) {m_pCamera->update(x);}
    // if (m_pOrbitCam) {m_pOrbitCam->update(x);}
}

void Program::render(int width, int height)
{
    mat4 mProj = m_pCamera->getProjMatrix(width, height);
    mat4 mView = m_pCamera->getViewMatrix();
    // mat4 mProj = m_pOrbitCam->getProjMatrix(width, height);
    // mat4 mView = m_pOrbitCam->getViewMatrix();

    city.render(mProj,mView,m_timer);
}

void Program::PressOne()
{
    city.reGenerate();
}

void Program::PressTwo()
{
    // city.deGenerate();
}