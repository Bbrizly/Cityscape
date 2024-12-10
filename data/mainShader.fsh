in vec4 v_color;      // Vertex color
in vec2 v_uv1;        // Texture coordinates
in float v_layer;     // Texture layer index
in vec3 v_normal;     // Vertex normal

out vec4 PixelColor;

uniform sampler2DArray u_texture;      // Texture array: layer 0 - building, layer 4 - window mask

uniform vec3 u_lightDir;               // Directional light direction (should be normalized)
uniform vec3 u_lightColor;             // Light color/intensity (e.g., white: vec3(1.0))
uniform vec3 u_ambient;                // Ambient light color/intensity (e.g., vec3(0.2))
uniform float u_windowLitFactor;       // Controls the probability of windows being lit

float windowMask = 0.0;          // Initialize window mask
float vcolorMask = 1.0;          // Mask for v_color influence (1.0 = full influence, 0.0 = no influence)

void main()
{
    vec4 baseColor = texture(u_texture, vec3(v_uv1, v_layer));
    vec3 finalColor = baseColor.rgb; // Start with the building texture color
    bool shouldLight = false;
    switch (int(v_layer))
    {
        case 0: case 1: case 2: //building 0,1,2 [+ 3] = red textures 3,4,5
            windowMask = texture(u_texture, vec3(v_uv1, v_layer + 3)).r; // Window mask from layer 4

            //Curtain logic
            vec2 scaledUV = floor(v_uv1 * ((v_layer + 5.0f) - (v_layer*2))); // * uvScaleMultiplier)); // Scale factor for window "cells

            float randVal = fract(sin(dot(scaledUV, vec2(12.9898, 78.233))) * 43758.5453);

            shouldLight = (windowMask > 0.5) && (randVal < u_windowLitFactor);

            vcolorMask = 1.0 - step(0.5, windowMask); // vcolorMask = 0.0 in window areas, 1.0 elsewhere
            break;
    }
    
    if (shouldLight)
    {
        vec3 windowColor = vec3(1.0, 1.0, 0.8); // Window light color (soft yellow)
        
        // Override the final color with the bright window light
        finalColor = mix(baseColor.rgb, windowColor, 0.4); // Stronger mix for bright windows
    }
    else
    {
        // Lightinng
        vec3 N = normalize(v_normal);  // Normalized surface normal
        vec3 L = normalize(u_lightDir); // Normalized light direction

        // Calculate diffuse lighting (Lambertian reflectance)
        float diff = max(dot(N, L), 0.0);

        // Combine ambient and diffuse lighting for the building
        vec3 lighting = u_ambient + (u_lightColor * diff);

        finalColor *= lighting;
    }

    finalColor = mix(finalColor, v_color.rgb, vcolorMask * 0.2); // Slightly influenced by v_color


    float finalAlpha = baseColor.a * v_color.a; // Combine alpha channels
    PixelColor = vec4(finalColor, finalAlpha);
}
