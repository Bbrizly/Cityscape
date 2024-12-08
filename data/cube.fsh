// Inputs from the vertex shader
in vec4 v_color;    // Vertex color
in vec2 v_uv1;      // Texture coordinates
in float v_layer;  
in vec3 v_normal;   // Vertex normal 

// Output to the framebuffer
out vec4 PixelColor;

// Uniforms
uniform sampler2DArray u_texture;

uniform vec3 u_lightDir;    // Directional light direction (should be normalized)
uniform vec3 u_lightColor;  // Light color/intensity (e.g., white: vec3(1.0))
uniform vec3 u_ambient;     // Ambient light color/intensity (e.g., vec3(0.2))

void main()
{
    vec3 N = normalize(v_normal);
    vec3 L = normalize(u_lightDir);

    float diff = max(dot(N, L), 0.0);

    // Sample the array texture using (u, v, layer)
    vec3 texCoords = vec3(v_uv1, v_layer);
    vec4 texColor = texture(u_texture, texCoords);

    vec3 lighting = u_ambient + (u_lightColor * diff);

    PixelColor = vec4(texColor.rgb * v_color.rgb * lighting, texColor.a * v_color.a);
    
}
