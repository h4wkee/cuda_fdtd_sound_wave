#pragma once

#include <mutex>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

struct GLFWwindow;

class Window
{
public:
	Window(unsigned width, unsigned height);
	~Window();

	bool update();
	bool isOpened(){return _opened;}

private:
	GLFWwindow * _window;

	bool _init(int width, int height);

	bool _opened = false;
	unsigned int _width;
	unsigned int _height;
};
