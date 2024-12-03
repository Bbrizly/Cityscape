
in vec4 v_color;
in vec2 v_uvcoord;
out vec4 PixelColor;

uniform sampler2D u_texture;

void main()
{
    vec4 texColor = texture(u_texture, v_uvcoord);
    PixelColor = texColor;// * v_color; 

    //PixelColor = vec4(v_uvcoord, 0.0, 1.0); // Red = U, Green = V
}