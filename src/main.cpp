#include <stdio.h>
#include <iostream>
#include <glm/glm.hpp>
#include "../wolf/wolf.h"
#include "../samplefw/SampleRunner.h"
#include "Program.h"

class Week2: public wolf::App
{
private:
    SampleRunner m_sampleRunner;
    Program* program;
public:
    Week2() : App("Citscape - Bassam")
    {
        program = new Program(this);
        m_sampleRunner.addSample(new Program(this));

        
    }

    ~Week2()
    {
    }

    void update(float dt) override
    {
        m_sampleRunner.update(dt);
        
        if(isKeyJustDown('1') || isKeyJustDown('r'))
        {
            program->PressOne();
        }
    }

    void render() override
    {
        m_sampleRunner.render(m_width, m_height);
    }
};

int main(int, char**) {
    Week2 week2;
    week2.run();
}