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

	GLFWmonitor * monitor = nullptr;
	if(!width) // full screen
	{
		monitor = glfwGetPrimaryMonitor();
		const GLFWvidmode * mode = glfwGetVideoMode(monitor);
		glfwWindowHint(GLFW_RED_BITS, mode->redBits);
		glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
		glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
		glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);
		_width = mode->width;
		_height = mode->height;
	}
	else
	{
		_width = width;
		_height = height;
	}

	glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
	glfwWindowHint(GLFW_DECORATED, GLFW_TRUE);
	glfwWindowHint(GLFW_MAXIMIZED, GLFW_FALSE);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_SAMPLES, 0);
	glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_TRUE);

	_window = glfwCreateWindow(_width, _height, "Cuda FDTD Sound Wave", monitor, nullptr);

	if (!_window)
	{
		std::cout << "Failed to open window" << std::endl;
		glfwTerminate();
		return false;
	}

	glfwMakeContextCurrent(_window);
	//glfwSwapInterval(true); // vsync enabled

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

	glfwMaximizeWindow(_window);

	glClear(GL_COLOR_BUFFER_BIT);
	glfwSwapBuffers(_window);

	// alpha blending
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

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
		glfwDestroyWindow(_window);
		glfwTerminate();
		return false;
	}
	return true;
}
