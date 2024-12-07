uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

in vec3 a_position;
in vec4 a_color;
in vec2 a_uvcoord;
in vec3 a_normal;

out vec4 v_color;
out vec2 v_uvcoord;
out vec3 v_normal;

void main()
{
    gl_Position = projection * view * world * vec4(a_position,1.0f);
	v_color = a_color;
    v_uvcoord = a_uvcoord;
    // v_normal = a_normal;
}