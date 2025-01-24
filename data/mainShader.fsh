in vec4 v_color;      // Vertex color
in vec2 v_uv1;        // Texture coordinates
in float v_layer;     // Texture layer index
in vec3 v_normal;     // Vertex normal
in vec3 v_worldPos;

out vec4 PixelColor;

uniform sampler2DArray u_texture;      // Texture array: layer 0 - building, layer 4 - window mask

uniform vec3 u_lightDir;               // Directional light direction (should be normalized)
uniform vec3 u_lightColor;             // Light color/intensity (e.g., white: vec3(1.0))
uniform vec3 u_ambient;                // Ambient light color/intensity (e.g., vec3(0.2))
uniform float u_windowLitFactor;       // Controls the probability of windows being lit

// Specular uniforms
uniform vec3 u_viewPos;           // Camera position in world space
uniform float u_specularStrength; // Specular highlight strength
uniform float u_shininess;        // Shininess factor

float windowMask = 0.0;          // Initialize window mask
float vcolorMask = 1.0;          // Mask for v_color influence (1.0 = full influence, 0.0 = no influence)

void main()
{
    vec4 baseColor = texture(u_texture, vec3(v_uv1, v_layer));
    vec3 finalColor = baseColor.rgb; // Start with the building texture color
    bool shouldLight = false;
    bool window = false;
    switch (int(v_layer))
    {
        case 0: case 1: case 2: //building 0,1,2 [+ 3] = red textures 3,4,5
            windowMask = texture(u_texture, vec3(v_uv1, v_layer + 3)).r; // Window mask from layer 4

            //Curtain logic
            vec2 scaledUV = floor(v_uv1 * ((v_layer + 5.0f) - (v_layer*2))); // Scale factor for window cells

            float randVal = fract(sin(dot(scaledUV, vec2(12.9898, 78.233))) * 43758.5453);

            window = windowMask > 0.5;
            shouldLight = window && (randVal < u_windowLitFactor);

            vcolorMask = 1.0 - step(0.5, windowMask); // vcolorMask = 0.0 in window areas, 1.0 elsewhere
            break;
    }

    vec3 N = normalize(v_normal);
    vec3 L = normalize(u_lightDir);
    vec3 V = normalize(u_viewPos - v_worldPos); // View direction

    // Diffuse + Ambient
    float diff = max(dot(N, L), 0.0);
    vec3 lighting = u_ambient + (u_lightColor * diff);
    
    if (shouldLight)
    {
        vec3 windowColor = vec3(1.0, 1.0, 0.8); // soft yellow
        
        // Override the final color with the bright window light
        finalColor = mix(baseColor.rgb, windowColor, 0.4); //mix for bright windows
    }
    else
    {
        // Lightinng
        vec3 N = normalize(v_normal);
        vec3 L = normalize(u_lightDir);

        float diff = max(dot(N, L), 0.0);

        vec3 lighting = u_ambient + (u_lightColor * diff);

        finalColor *= lighting;
    }
    
    if (window) {
        vec3 R = reflect(-L, N);
        float spec = pow(max(dot(V, R), 0.0), u_shininess);
        vec3 specular = u_lightColor * u_specularStrength * spec;
        finalColor += specular;
    }

    finalColor = mix(finalColor, v_color.rgb, vcolorMask * 0.2); // Slightly influenced by v_color


    float finalAlpha = baseColor.a * v_color.a;
    PixelColor = vec4(finalColor, finalAlpha);
}
