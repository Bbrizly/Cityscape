#pragma once
#include "../wolf/wolf.h"
#include "../samplefw/Sample.h"
#include "../samplefw/OrbitCamera.h"
#include "../samplefw/FirstPersonCamera.h"

class Program: public Sample
{
private:
    GLuint cubeShader = 0;
    GLuint otherShader = 0;
    wolf::Program* m_pProgram = 0;
    GLuint m_vbo = 0;
    GLuint m_vao = 0; 
    float m_timer = 0;
    FirstPersonCamera* m_pCamera = nullptr;
public:
    Program(wolf::App* pApp) : Sample(pApp,"Plane") {}
    ~Program();
    
    void PressOne();
    void PressTwo();

    void init() override;
    void update(float dt) override;
    void render(int width, int height) override;
    
};
