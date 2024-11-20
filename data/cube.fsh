in vec4 v_color;
out vec4 PixelColor;

//del
uniform vec3 objectColor;

void main()
{
    //PixelColor = vec4(objectColor, 1.0);
    PixelColor = v_color;
}