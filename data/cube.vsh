uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;

in vec3 a_position;
in vec4 a_color;
in vec2 a_uv1;
in float a_layer; // New: layer index
in vec3 a_normal;

out vec4 v_color;
out vec2 v_uv1;
out float v_layer; // Pass to fragment shader
out vec3 v_normal;

void main()
{
    gl_Position = projection * view * world * vec4(a_position,1.0f);
	v_color = a_color;
    v_uv1 = a_uv1;
    // v_normal = vec4(a_normal,1.0f);
    
    mat3 normalMatrix = mat3(transpose(inverse(world)));
    vec3 transformedNormal = normalize(normalMatrix * a_normal);
    // v_normal = transformedNormal;
    v_normal = a_normal;
    
    v_layer = a_layer;

    // v_normal = mat3(transpose(inverse(world))) * a_normal; // Correct normal transformation
}