/*
uniform mat4 view;
uniform mat4 projection;
uniform mat4 world;

in vec3 aPos;
in vec4 aColor;
in vec4 instanceData;
in vec4 a_color;

out vec4 v_color;

void main() {
    vec3 instancePos = instanceData.xyz;
    float scale = instanceData.w;

    // Scale the building height
    vec3 scaledPos = vec3(aPos.x, aPos.y * scale, aPos.z);

    vec4 worldPos = vec4(scaledPos + instancePos, 1.0);

    gl_Position = projection * view * worldPos;
    v_color = a_color;
}
*/

uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

in vec4 a_position;

in vec4 a_color;

out vec4 v_color;

void main()
{
    gl_Position = projection * view * world * a_position;
	v_color = a_color;
}