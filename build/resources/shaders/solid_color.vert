layout (location = 0) in vec3 position;
layout (location = 1) in vec3 color;

out vec3 Color;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

void main()
{
    gl_Position = projection * view * model * vec4(position, 1.0);
    Color = color;
    //gl_Position = vec4(position, 1.0);
}
