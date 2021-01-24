#define _USE_MATH_DEFINES
#include <iostream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <chrono>

#include <Window.h>
#include <Shader.h>
#include <Renderable.h>
#include <AcousticFDTD.h>

const float windowWidth = 800.f;
const float windowHeight = 600.f;

int main(int argc, char * argv[])
{
	srand(time(NULL));
	/////////////// FDTD INITIALIZATION
	float surfaceScale = 5.f;
	glm::vec2 surfaceSize = {2 * surfaceScale, 2 * surfaceScale};
	float resolution = 0.1f;
	bool oglContext = false;
	for(unsigned int i = 0; i < argc; ++i)
	{
		if(argc > i + 1)
		{
			if(!strcmp(argv[i], "-r"))
			{
				resolution = 10.f/atof(argv[2]);
			}
		}
		if(!strcmp(argv[i], "-ogl"))
		{
			oglContext = true;
		}
	}
	glm::ivec2 gridSize = {surfaceSize.x / resolution, surfaceSize.y / resolution};
	glm::vec3 * gridColors = new glm::vec3[gridSize.x * gridSize.y];
	AcousticFDTD fdtd(gridColors, gridSize);
	///////////////

	if(oglContext)
	{
		Window window(windowWidth, windowHeight);
		if(!window.isOpened())
		{
			std::cout << "Couldn't create a window!" << std::endl;
			return 1;
		}

		/////////////// OPENGL INITIALIZATION
		float pointSize = resolution * 100.f;
		float renderPadding = 0.1f;
		glClearColor(0.33f, 0.33f, 0.33f, 1.f);

		Shader shader("solid_color");
		//glm::mat4 projection = glm::ortho(0.0f, windowWidth, 0.0f, windowHeight, 0.1f, 100.f);
		glm::mat4 projection = glm::perspective(glm::radians(45.f), windowWidth / windowHeight, 0.1f, 100.f);
		glm::mat4 view = view = glm::lookAt(glm::vec3(0.0f, -10.0f, -10.0f),glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
		glm::mat4 model = glm::mat4(1.f);

		glm::vec3 surfaceColor = {0.7f, 0.7f, 0.7f};
		glm::vec3 pointsColor = {1.f, 1.f, 1.f};

		std::vector<glm::vec3> surfacePositions = {
				{-1.f * surfaceScale,1.f * surfaceScale,1.f},
				{1.f * surfaceScale, 1.f * surfaceScale,1.f},
				{1.f * surfaceScale, -1.f * surfaceScale, 1.f},
				{1.f * surfaceScale, -1.f * surfaceScale, 1.f},
				{-1.f * surfaceScale, -1.f * surfaceScale, 1.f},
				{-1.f * surfaceScale, 1.f * surfaceScale, 1.f}
		};
		Renderable surface(surfacePositions, GL_TRIANGLES);

		Renderable * points = new Renderable[gridSize.x * gridSize.y];

		std::vector<glm::vec3> pointPosition;
		pointPosition.resize(1);
		glm::vec2 base = {-(gridSize.x - resolution) * resolution / 2.f, -(gridSize.y - resolution) * resolution / 2.f};
		for(unsigned int i = 0; i < gridSize.x; ++i)
		{
			for(unsigned j = 0; j < gridSize.y; ++j)
			{
				pointPosition[0] = {i * resolution + base.x, j * resolution + base.y, 1.f};
				points[i * gridSize.y + j].init(pointPosition, GL_POINTS, pointSize);
			}
		}

		//	std::vector<glm::vec3> pointsPositions = {
//			{-0.5f,0.5f,1.f},
//			{0.5f, 0.5f,1.f},
//			{0.f, -0.5f, 1.f}
//	};

		//Renderable points(pointsPositions, GL_POINTS);

		/////////////////

		auto appStart = std::chrono::high_resolution_clock::now();

		while(window.update())
		{
			auto appTime = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(appTime - appStart).count();
			pointsColor = glm::vec3{sin(duration * 0.01f), 0.5f, 1.f};

			shader.setUniform4m("projection", glm::value_ptr(projection));
			shader.setUniform4m("view", glm::value_ptr(view));
			shader.setUniform4m("model", glm::value_ptr(model));

			shader.setUniform3("color", surfaceColor.x, surfaceColor.y, surfaceColor.z);
			surface.draw();

//		shader.setUniform3("color", pointsColor.x, pointsColor.y, pointsColor.z);
//		points.draw();

			fdtd.draw();
			for(unsigned int i = 0; i < gridSize.x; ++i)
			{
				for(unsigned j = 0; j < gridSize.y; ++j)
				{
					glm::vec3 pointColor = gridColors[i * gridSize.y + j];
					//std::cout << pointColor.x << ',' << pointColor.y << ',' << pointColor.z << std::endl;
					shader.setUniform3("color", pointColor.x, 1.f, 1.f);
					points[i * gridSize.y + j].draw();
				}
			}
		}

		delete[] gridColors;
	}
	else
	{
		AcousticFDTD * fdtd = new AcousticFDTD(nullptr, gridSize);
		while(true)
		{
			fdtd->draw();
		}
	}
	return 0;
}
