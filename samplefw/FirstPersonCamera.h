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

    float m_normalSpeed = 2.0f;
    float m_SprintSpeed = 5.0f;

    glm::vec2 m_lastMousePos; // Store the last mouse position
    float m_movementSpeed;    // Movement speed
    float m_mouseSensitivity; // Mouse sensitivity
};

#endif // FIRSTPERSONCAMERA_H


/*#pragma once

#include "../wolf/wolf.h"

class FirstPersonCamera
{
public:
    FirstPersonCamera(wolf::App* pApp);
    virtual ~FirstPersonCamera();

    void update(float dt);
    glm::mat4 getViewMatrix();
    glm::mat4 getProjMatrix(int width, int height);
    glm::vec3 getViewDirection() const;
    glm::vec3 getViewPosition() const;

    // void focusOn(const glm::vec3& min, const glm::vec3& max);

private:
    // void _rotate(const glm::vec2& mouseMovement);
    // glm::vec3 _getCameraUp();
    // glm::vec3 _getCameraSide();
    // void _pan(const glm::vec2& mouseMovement);
    // float _calculateRequiredDistance();

    // float m_rotX                = 0.0f;
    // float m_rotY                = 0.0f;
    // float m_distance            = 100.0f;
    glm::vec3 m_offset          = glm::vec3(0.0f,0.0f,0.0f);
    glm::vec3 m_position        = glm::vec3(0.0f,0.0f,0.0f);
    glm::vec3 m_target          = glm::vec3(0.0f,0.0f,0.0f);
    glm::vec3 m_focusMin        = glm::vec3(0.0f,0.0f,0.0f);
    glm::vec3 m_focusMax        = glm::vec3(0.0f,0.0f,0.0f);
    float m_fov                 = glm::radians(190.0f);
    float m_near                = 0.1f;
    float m_far                 = 1000.0f;

    float m_pitch; // Rotation around X-axis
    float m_yaw;   // Rotation around Y-axis

    glm::vec2 m_lastMousePos    = glm::vec2(0.0f,0.0f);
    wolf::App* m_pApp           = nullptr;
};*/