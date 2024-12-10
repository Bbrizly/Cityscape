#ifndef FIRSTPERSONCAMERA_H
#define FIRSTPERSONCAMERA_H

#include "../wolf/wolf.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class FirstPersonCamera
{
public:
    FirstPersonCamera(wolf::App* pApp);
    ~FirstPersonCamera();

    void update(float dt);
    glm::mat4 getViewMatrix() const;
    glm::mat4 getProjMatrix(int width, int height) const;

    glm::vec3 getPosition() const { return m_position; }
    glm::vec3 getDirection() const { return m_direction; }

private:
    void _updateOrientation(const glm::vec2& mouseMovement);

    wolf::App* m_pApp;

    glm::vec3 m_position = glm::vec3(0.0f,50.0f,0.0f);  // Camera position
    glm::vec3 m_direction; // Forward direction
    glm::vec3 m_up;        // Up vector
    float m_yaw;           // Rotation around the Y-axis
    float m_pitch;         // Rotation around the X-axis
    bool m_invertY;

    float m_normalSpeed = 20.0f;
    float m_SprintSpeed = 50.0f;

    glm::vec2 m_lastMousePos; // Store the last mouse position
    float m_movementSpeed;    // Movement speed
    float m_mouseSensitivity; // Mouse sensitivity
};

#endif