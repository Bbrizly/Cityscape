#include "FirstPersonCamera.h"

FirstPersonCamera::FirstPersonCamera(wolf::App* pApp)
    : m_pApp(pApp),
      m_position(0.0f, 0.0f, 5.0f),
      m_direction(0.0f, 0.0f, -1.0f),
      m_up(0.0f, 1.0f, 0.0f),
      m_yaw(-90.0f),
      m_pitch(0.0f),
      m_movementSpeed(5.0f),
      m_mouseSensitivity(0.1f)
{
    m_lastMousePos = m_pApp->getMousePos();
}

FirstPersonCamera::~FirstPersonCamera() {}

void FirstPersonCamera::update(float dt)
{
    float speedMultiplier = m_pApp->isKeyDown(GLFW_KEY_LEFT_SHIFT) ? 2.5f : 1.0f;
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

    // Handle mouse movement for looking around
    glm::vec2 currentMousePos = m_pApp->getMousePos();
    glm::vec2 mouseMovement = currentMousePos - m_lastMousePos;
    m_lastMousePos = currentMousePos;

    _updateOrientation(mouseMovement);
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



/*#include "FirstPersonCamera.h"

const float MOVEMENT_SPEED = 2.0f;
const float MOUSE_SENSITIVITY = 3.0f;

FirstPersonCamera::FirstPersonCamera(wolf::App* pApp)
    : m_pApp(pApp), m_position(0.0f, 0.0f, 5.0f), // Starting position
      m_pitch(0.0f), m_yaw(-90.0f)
{
    m_lastMousePos = m_pApp->getMousePos();
}

FirstPersonCamera::~FirstPersonCamera() {}

void FirstPersonCamera::update(float dt)
{
    // Handle keyboard input for movement
    glm::vec3 front;
    front.x = cos(glm::radians(m_pitch)) * cos(glm::radians(m_yaw));
    front.y = sin(glm::radians(m_pitch));
    front.z = cos(glm::radians(m_pitch)) * sin(glm::radians(m_yaw));
    front = glm::normalize(front);

    glm::vec3 right = glm::normalize(glm::cross(front, glm::vec3(0.0f, 1.0f, 0.0f)));
    glm::vec3 up = glm::normalize(glm::cross(right, front));

    // Movement: W, A, S, D
    if (m_pApp->isKeyDown('W'))
    {
        m_position += front * MOVEMENT_SPEED * dt;
    }
    if (m_pApp->isKeyDown('S'))
    {
        m_position -= front * MOVEMENT_SPEED * dt;
    }
    if (m_pApp->isKeyDown('A'))
    {
        m_position -= right * MOVEMENT_SPEED * dt;
    }
    if (m_pApp->isKeyDown('D'))
    {
        m_position += right * MOVEMENT_SPEED * dt;
    }

    // Handle mouse movement for camera rotation when RMB is held down
    if (m_pApp->isRMBDown())
    {
        glm::vec2 mousePos = m_pApp->getMousePos();
        glm::vec2 mouseMovement = mousePos - m_lastMousePos;

        m_yaw += mouseMovement.x * MOUSE_SENSITIVITY;
        m_pitch -= mouseMovement.y * MOUSE_SENSITIVITY;

        // Constrain the pitch angle to prevent screen flip
        if (m_pitch > 89.0f)
            m_pitch = 89.0f;
        if (m_pitch < -89.0f)
            m_pitch = -89.0f;
    }

    m_lastMousePos = m_pApp->getMousePos();
}

glm::mat4 FirstPersonCamera::getViewMatrix()
{
    glm::vec3 front;
    front.x = cos(glm::radians(m_pitch)) * cos(glm::radians(m_yaw));
    front.y = sin(glm::radians(m_pitch));
    front.z = cos(glm::radians(m_pitch)) * sin(glm::radians(m_yaw));
    front = glm::normalize(front);

    return glm::lookAt(m_position, m_position + front, glm::vec3(0.0f, 1.0f, 0.0f));
}

glm::mat4 FirstPersonCamera::getProjMatrix(int width, int height)
{
    return glm::perspective(glm::radians(m_fov), (float)width / (float)height, m_near, m_far);}
*/
/* void FirstPersonCamera::focusOn(const glm::vec3& min, const glm::vec3& max)
// {
//     m_focusMin = min;
//     m_focusMax = max;
//     m_offset = glm::vec3(0.0f,0.0f,0.0f);
//     m_rotX = -MATH_PI / 4.0f;
//     m_rotY = MATH_PI / 4.0f;

//     m_target = min + ((max - min) * 0.5f);

//     m_distance = _calculateRequiredDistance();
// }

// void FirstPersonCamera::_rotate(const glm::vec2& mouseMovement)
// {
//     m_rotX -= mouseMovement.y * 0.003f;
//     m_rotY -= mouseMovement.x * 0.003f;
// }

// glm::vec3 FirstPersonCamera::_getCameraSide()
// {
//     glm::vec3 dir = m_target - m_position;
//     glm::vec3 side = glm::cross(dir, glm::vec3(0.0f,1.0f,0.0f));
//     return glm::normalize(side);
// }

// glm::vec3 FirstPersonCamera::_getCameraUp()
// {
//     glm::vec3 dir = m_target - m_position;
//     glm::vec3 v = _getCameraSide();
//     v = glm::cross(dir, v);
//     return glm::normalize(v);
// }

// void FirstPersonCamera::_pan(const glm::vec2& mouseMovement)
// {
//     glm::vec3 side = _getCameraSide();
//     glm::vec3 up = _getCameraUp();

//     side = side * -mouseMovement[0] * 0.007f * (m_distance / 5.0f);
//     up = up * -mouseMovement[1] * 0.007f * (m_distance / 5.0f);

//     m_offset += side;
//     m_offset += up;
// }    

// float FirstPersonCamera::_calculateRequiredDistance() 
// {
//     glm::vec3 min = m_focusMin;
//     glm::vec3 max = m_focusMax;
//     glm::vec3 center = min + ((max - min) * 0.5f);
//     float r = wolf::max(glm::distance(center,min), glm::distance(center,max));

//     return (r * 2.0f) / tan(m_fov / 1.5f);
// }

glm::vec3 FirstPersonCamera::getViewDirection() const
{
    glm::vec3 front;
    front.x = cos(glm::radians(m_pitch)) * cos(glm::radians(m_yaw));
    front.y = sin(glm::radians(m_pitch));
    front.z = cos(glm::radians(m_pitch)) * sin(glm::radians(m_yaw));
    return glm::normalize(front);
}

glm::vec3 FirstPersonCamera::getViewPosition() const
{
    return m_position;
}

*/
/*
#include "FirstPersonCamera.h"

FirstPersonCamera::FirstPersonCamera(wolf::App* pApp)
    : m_pApp(pApp)
{
    m_lastMousePos = m_pApp->getMousePos();
}

FirstPersonCamera::~FirstPersonCamera()
{

}

void FirstPersonCamera::update(float dt)
{
    glm::vec2 mousePos = m_pApp->getMousePos();

    if(m_pApp->isLMBDown())
    {
        glm::vec2 mouseMovement = mousePos - m_lastMousePos;
        _rotate(mouseMovement);
    }
    else if(m_pApp->isMMBDown())
    {
        glm::vec2 mouseMovement = mousePos - m_lastMousePos;
        _pan(mouseMovement);
    }

    glm::vec2 mouseScroll = m_pApp->getMouseScroll();

    if(mouseScroll.y > 0) {
        m_distance -= (m_distance / 5.0f);
    } else if(mouseScroll.y < 0) {
        m_distance += (m_distance / 5.0f);
    }

    m_distance = wolf::max(10.0f, m_distance);

    m_far = wolf::max(150.0f, m_distance * 2.0f);
    m_near = m_distance / 10.0f;

    if(m_pApp->isKeyDown('f'))
    {
        focusOn(m_focusMin,m_focusMax);
    }

    m_lastMousePos = mousePos;
}

glm::mat4 FirstPersonCamera::getViewMatrix()
{
    glm::mat4 m(1.0f);
    m = glm::rotate(m_rotY, glm::vec3(0.0f,1.0f,0.0f));
    m = m * glm::rotate(m_rotX, glm::vec3(1.0f,0.0f,0.0f));

    m_position = m * glm::vec4(0.0f,0.0f,m_distance,1.0f);

    glm::vec3 up = m * glm::vec4(0.0f,1.0f,0.0f,1.0f);
    glm::vec3 pos = m_position + m_offset;
    glm::vec3 target = m_target + m_offset;

    return glm::lookAt(pos, target, up);
}

glm::mat4 FirstPersonCamera::getProjMatrix(int width, int height)
{
	return glm::perspective(m_fov, (float)width / (float)height, m_near, m_far);
}

void FirstPersonCamera::focusOn(const glm::vec3& min, const glm::vec3& max)
{
    m_focusMin = min;
    m_focusMax = max;
    m_offset = glm::vec3(0.0f,0.0f,0.0f);
    m_rotX = -MATH_PI / 4.0f;
    m_rotY = MATH_PI / 4.0f;

    m_target = min + ((max - min) * 0.5f);

    m_distance = _calculateRequiredDistance();
}

void FirstPersonCamera::_rotate(const glm::vec2& mouseMovement)
{
    m_rotX -= mouseMovement.y * 0.003f;
    m_rotY -= mouseMovement.x * 0.003f;
}

glm::vec3 FirstPersonCamera::_getCameraSide()
{
    glm::vec3 dir = m_target - m_position;
    glm::vec3 side = glm::cross(dir, glm::vec3(0.0f,1.0f,0.0f));
    return glm::normalize(side);
}

glm::vec3 FirstPersonCamera::_getCameraUp()
{
    glm::vec3 dir = m_target - m_position;
    glm::vec3 v = _getCameraSide();
    v = glm::cross(dir, v);
    return glm::normalize(v);
}

void FirstPersonCamera::_pan(const glm::vec2& mouseMovement)
{
    glm::vec3 side = _getCameraSide();
    glm::vec3 up = _getCameraUp();

    side = side * -mouseMovement[0] * 0.007f * (m_distance / 5.0f);
    up = up * -mouseMovement[1] * 0.007f * (m_distance / 5.0f);

    m_offset += side;
    m_offset += up;
}    

float FirstPersonCamera::_calculateRequiredDistance() 
{
    glm::vec3 min = m_focusMin;
    glm::vec3 max = m_focusMax;
    glm::vec3 center = min + ((max - min) * 0.5f);
    float r = wolf::max(glm::distance(center,min), glm::distance(center,max));

    return (r * 2.0f) / tan(m_fov / 1.5f);
}

glm::vec3 FirstPersonCamera::getViewDirection() const
{
    return glm::normalize(m_target - m_position);
}

glm::vec3 FirstPersonCamera::getViewPosition() const
{
    return m_position;
}
*/