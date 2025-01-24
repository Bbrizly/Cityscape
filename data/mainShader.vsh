uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

in vec3 a_position;
in vec4 a_color;
in vec2 a_uv1;
in float a_uv2;
in vec3 a_normal;

out vec4 v_color;
out vec2 v_uv1;
out float v_layer;
out vec3 v_normal;
out vec3 v_worldPos;

void main()
{
    vec4 worldPos = world * vec4(a_position, 1.0);
    gl_Position = projection * view * worldPos;
    v_color = a_color;
    v_uv1 = a_uv1;

    mat3 normalMatrix = mat3(transpose(inverse(world)));
    vec3 transformedNormal = normalize(normalMatrix * a_normal);
    v_normal = transformedNormal;

    v_layer = a_uv2;
    v_worldPos = worldPos.xyz;
}