layout (location = 0) out vec4 outColor;

in vec3 Color;

//uniform vec3 color;

void main()
{
	//RESOLVE:
	outColor = vec4(Color, 1.0);
}
