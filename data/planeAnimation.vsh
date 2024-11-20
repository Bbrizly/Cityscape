uniform mat4 projection;
uniform mat4 view;
uniform mat4 world;
uniform float u_time;

uniform float speed;        // 5.0f
uniform float frequency;    // 5.0f
uniform float amplitude;    // 0.5f
uniform bool chooseColor;   // true = normal blue water waves, false = rainbow

in vec3 a_position;
in vec4 a_color;

out vec4 v_color;

vec4 getRainbowColor(float t) {
    // Clamp t around 1.0 to ensure it stays in the [0, 1] range
    t = mod(t, 1.0);

    if (t < 1.0 / 6.0) {
        return vec4(1.0, t * 6.0, 0.0, 1.0); // Red to Yellow
    } else if (t < 2.0 / 6.0) {
        return vec4(1.0 - (t - 1.0 / 6.0) * 6.0, 1.0, 0.0, 1.0); // Yellow to Green
    } else if (t < 3.0 / 6.0) {
        return vec4(0.0, 1.0, (t - 2.0 / 6.0) * 6.0, 1.0); // Green to Cyan
    } else if (t < 4.0 / 6.0) {
        return vec4(0.0, 1.0 - (t - 3.0 / 6.0) * 6.0, 1.0, 1.0); // Cyan to Blue
    } else if (t < 5.0 / 6.0) {
        return vec4((t - 4.0 / 6.0) * 6.0, 0.0, 1.0, 1.0); // Blue to Magenta
    } else {
        return vec4(1.0, 0.0, 1.0 - (t - 5.0 / 6.0) * 6.0, 1.0); // Magenta to Red
    }
}

void main() {
    vec4 pos = vec4(a_position, 1.0);

    // Animate the y-coordinate using a sine wave function based on both position and time
    float wave = sin(pos.x * frequency + u_time * speed) * 0.5;
    pos.z += wave * amplitude;  // Adjust wave amplitude

    gl_Position = projection * view * world * pos;

    float colorModifier = (wave + 0.1) / 1.0; // Normalize wave to be in range [0, 1]
    
    vec4 baseColor = vec4(0.0, 0.0, 0.5, 1.0);  // dark blue color
    vec4 topColor = vec4(0.4, 0.4, 1.0, 1.0);   // Light blue at the peak of the wave
    vec4 animatedColor = mix(baseColor, topColor, colorModifier);

    if (chooseColor) {
        v_color = animatedColor;      
    } else {
        v_color = getRainbowColor(colorModifier);
    }
}
