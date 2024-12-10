#include "FirstPersonCamera.h"
#include <glm/gtx/norm.hpp> // for length2 if needed

FirstPersonCamera::FirstPersonCamera(wolf::App* pApp)
    : m_pApp(pApp),
      m_position(0.0f, 50.0f, 5.0f),
      m_direction(0.0f, 0.0f, -1.0f),
      m_up(0.0f, 1.0f, 0.0f),
      m_yaw(-90.0f),
      m_pitch(0.0f),
      m_movementSpeed(10.0f),
      m_mouseSensitivity(0.1f),
      m_invertY(false), // New: invert Y toggle
      m_normalSpeed(10.0f),
      m_SprintSpeed(20.0f)
{
    m_lastMousePos = m_pApp->getMousePos();
}

FirstPersonCamera::~FirstPersonCamera() {}

void FirstPersonCamera::update(float dt)
{
    float speedMultiplier = m_pApp->isKeyDown(GLFW_KEY_LEFT_SHIFT) ? m_SprintSpeed : m_normalSpeed;
    float adjustedSpeed = m_movementSpeed * speedMultiplier;

    // Compute basis vectors
    glm::vec3 right = glm::normalize(glm::cross(m_direction, m_up));
    glm::vec3 forward = glm::normalize(m_direction);
    glm::vec3 movement(0.0f);

    // Check movement keys
    if (m_pApp->isKeyDown('W')) {
        movement += forward;
    }
    if (m_pApp->isKeyDown('S')) {
        movement -= forward;
    }
    if (m_pApp->isKeyDown('A')) {
        movement -= right;
    }
    if (m_pApp->isKeyDown('D')) {
        movement += right;
    }

//space for up and lceft ctrl for down
    if (m_pApp->isKeyDown(GLFW_KEY_SPACE)) {
        movement += m_up;
    }
    if (m_pApp->isKeyDown(GLFW_KEY_LEFT_CONTROL)) {
        movement -= m_up;
    }

    // // Normalize movement so diagonals are not faster
    // float length = glm::length(movement);
    // if (glm::length(movement) > 0.0001f) {
    //     movement = glm::normalize(movement) * (adjustedSpeed * dt);
    // }
    // m_position += movement;
    movement *= adjustedSpeed * dt;

    // Update camera position
    m_position += movement;

    //Press Y to invert mouse y axis
    if (m_pApp->isKeyJustDown('Y')) {
        m_invertY = !m_invertY;
    }

    if (m_pApp->isRMBDown()) {
        glm::vec2 currentMousePos = m_pApp->getMousePos();
        glm::vec2 mouseMovement = currentMousePos - m_lastMousePos;
        _updateOrientation(mouseMovement);
        m_lastMousePos = currentMousePos; // Update only if RMB is held
    } else {
        m_lastMousePos = m_pApp->getMousePos(); 
    }
}

void FirstPersonCamera::_updateOrientation(const glm::vec2& mouseMovement)
{
    m_yaw += mouseMovement.x * m_mouseSensitivity;
    // Invert Y if m_invertY
    float yMovement = m_invertY ? mouseMovement.y : -mouseMovement.y;
    m_pitch += yMovement * m_mouseSensitivity; 

    // Clamp pitch
    m_pitch = glm::clamp(m_pitch, -89.0f, 89.0f);

    // Update direction
    glm::vec3 direction;
    direction.x = cos(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
    direction.y = sin(glm::radians(m_pitch));
    direction.z = sin(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
    m_direction = glm::normalize(direction);
}

glm::mat4 FirstPersonCamera::getViewMatrix() const
{
    return glm::lookAt(m_position, m_position + m_direction, m_up);
}

glm::mat4 FirstPersonCamera::getProjMatrix(int width, int height) const
{
    float fov = 60.0f;
    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    float nearPlane = 0.1f;
    float farPlane = 5000.0f;
    return glm::perspective(glm::radians(fov), aspectRatio, nearPlane, farPlane);
}


/*

#include "FirstPersonCamera.h"

FirstPersonCamera::FirstPersonCamera(wolf::App* pApp)
    : m_pApp(pApp),
      m_position(0.0f, 50.0f, 5.0f),
      m_direction(0.0f, 0.0f, -1.0f),
      m_up(0.0f, 1.0f, 0.0f),
      m_yaw(-90.0f),
      m_pitch(0.0f),
      m_movementSpeed(10.0f),
      m_mouseSensitivity(0.1f)
{
    m_lastMousePos = m_pApp->getMousePos();
}

FirstPersonCamera::~FirstPersonCamera() {}

void FirstPersonCamera::update(float dt)
{
    float speedMultiplier = m_pApp->isKeyDown(GLFW_KEY_LEFT_SHIFT) ? m_SprintSpeed : m_normalSpeed;
    float adjustedSpeed = m_movementSpeed * speedMultiplier;
    glm::vec3 right = glm::normalize(glm::cross(m_direction, m_up));
    glm::vec3 forward = glm::normalize(m_direction);

    if (m_pApp->isKeyDown('W')) {
        m_position += forward * adjustedSpeed * dt;
    }
    if (m_pApp->isKeyDown('S')) {
        m_position -= forward * adjustedSpeed * dt;
    }
    if (m_pApp->isKeyDown('A')) {
        m_position -= right * adjustedSpeed * dt;
    }
    if (m_pApp->isKeyDown('D')) {
        m_position += right * adjustedSpeed * dt;
    }
    if (m_pApp->isKeyDown(GLFW_KEY_SPACE)) {
        m_position += m_up * adjustedSpeed * dt;
    }
    if (m_pApp->isKeyDown(GLFW_KEY_LEFT_CONTROL)) {
        m_position -= m_up * adjustedSpeed * dt;
    }
    
    if (m_pApp->isRMBDown()) {
        glm::vec2 currentMousePos = m_pApp->getMousePos();
        glm::vec2 mouseMovement = currentMousePos - m_lastMousePos;
        _updateOrientation(mouseMovement);
        m_lastMousePos = currentMousePos; // Update the last mouse position only if RMB is held
    } else {
        m_lastMousePos = m_pApp->getMousePos(); // Track mouse position for smooth transition
    }
}

void FirstPersonCamera::_updateOrientation(const glm::vec2& mouseMovement)
{
    m_yaw += mouseMovement.x * m_mouseSensitivity;
    m_pitch -= mouseMovement.y * m_mouseSensitivity;

    // Clamp the pitch to avoid flipping the camera
    m_pitch = glm::clamp(m_pitch, -89.0f, 89.0f);

    // Update the direction vector based on yaw and pitch
    glm::vec3 direction;
    direction.x = cos(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
    direction.y = sin(glm::radians(m_pitch));
    direction.z = sin(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
    m_direction = glm::normalize(direction);
}

glm::mat4 FirstPersonCamera::getViewMatrix() const
{
    return glm::lookAt(m_position, m_position + m_direction, m_up);
}

glm::mat4 FirstPersonCamera::getProjMatrix(int width, int height) const
{
    float fov = 45.0f; // Field of view
    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    float nearPlane = 0.1f;
    float farPlane = 500.0f;
    return glm::perspective(glm::radians(fov), aspectRatio, nearPlane, farPlane);
}
*/