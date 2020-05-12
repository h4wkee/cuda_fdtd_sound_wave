#include <Window.h>

#include <iostream>

Window::Window(unsigned int width, unsigned int height)
{
	_init(width, height);
}

Window::~Window()
{

}

void escCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

void framebufferSizeCallback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void errorCallback(int error, const char* description)
{
	fprintf(stderr, "Error: %s\n", description);
}

bool Window::_init(int width, int height)
{
	if (!glfwInit())
	{
		return false;
	}

	glfwSetErrorCallback(errorCallback);

	_width = width;
	_height = height;

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

	_window = glfwCreateWindow(_width, _height, "Cuda FDTD Sound Wave", nullptr, nullptr);

	if (!_window)
	{
		std::cout << "Failed to open window" << std::endl;
		glfwTerminate();
		return false;
	}

	glfwMakeContextCurrent(_window);

	if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
	{
		std::cout << "Failed to open window: cannot initialize GLAD" << std::endl;
		return false;
	}

	if (!GLAD_GL_VERSION_3_3)
	{
		std::cout << "Failed to open window: OpenGL 3.3 or later is required, but " << GLVersion.major << "."
		          << GLVersion.minor << " was detected" << std::endl;
		return false;
	}

	std::cout << "OpenGL context version: " << GLVersion.major << "." << GLVersion.minor << std::endl;

	// store pointer to Window in GLFW window for later use
	glfwSetWindowUserPointer(_window, this);

	glClear(GL_COLOR_BUFFER_BIT);
	glfwSwapBuffers(_window);

	glfwSetKeyCallback(_window, escCallback);
	glfwSetFramebufferSizeCallback(_window, framebufferSizeCallback);
	_opened = true;

	return true;
}

bool Window::update()
{
	if (!glfwWindowShouldClose(_window))
	{
		glfwSwapBuffers(_window);
		glClear(GL_COLOR_BUFFER_BIT);
		glfwPollEvents();
	}
	else
	{
		return false;
	}
	return true;
}

void Window::close()
{
	if(_window)
	{
		glfwDestroyWindow(_window);
		glfwTerminate();
	}
}
