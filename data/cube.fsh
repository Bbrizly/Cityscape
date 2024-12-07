
in vec4 v_color;
in vec2 v_uvcoord;
in vec3 v_normal;
out vec4 PixelColor;

uniform sampler2D u_texture;

void main()
{
    //vec3 normalizedNormal = normalize(v_normal);
    //vec3 color = normalizedNormal * 0.5 + 0.5;
    //PixelColor = vec4(color, 1.0);

    vec4 texColor = texture(u_texture, v_uvcoord);
    PixelColor = texColor * v_color; 
    //PixelColor = v_color; 

    //PixelColor = vec4(v_uvcoord, 0.0, 1.0); // Red = U, Green = V
}
