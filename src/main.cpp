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
	Window window(windowWidth, windowHeight);
	if(!window.isOpened())
	{
		std::cout << "Couldn't create a window!" << std::endl;
		return 1;
	}
	/////////////// FDTD INITIALIZATION
	float surfaceScale = 5.f;
	glm::vec2 surfaceSize = {2 * surfaceScale, 2 * surfaceScale};
	float resolution = 0.1f;
	unsigned int dataPerThread = 1;
	unsigned int blockSize = 32; // 32 x 32
	for(unsigned int i = 1; i < argc; ++i)
	{
		if(argc > i + 1)
		{
			if (!strcmp(argv[i], "-r"))
			{
				resolution = atof(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-dpt"))
			{
				dataPerThread = atoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-bs"))
			{
				blockSize = atoi(argv[i + 1]);
			}
		}
	}
	glm::ivec2 gridSize = {surfaceSize.x / resolution, surfaceSize.y / resolution};
	///////////////

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

	std::vector<Vertex> surfaceVertices = {
			{
				{-1.f * surfaceScale, 1.f * surfaceScale, 1.f },
				{surfaceColor.x, surfaceColor.y, surfaceColor.z }
			},
			{
				{1.f * surfaceScale, 1.f * surfaceScale,1.f },
				{surfaceColor.x, surfaceColor.y, surfaceColor.z }
			},
			{
				{1.f * surfaceScale, -1.f * surfaceScale, 1.f },
				{surfaceColor.x, surfaceColor.y, surfaceColor.z }
			},
			{
				{1.f * surfaceScale, -1.f * surfaceScale, 1.f },
				{surfaceColor.x, surfaceColor.y, surfaceColor.z }
			},
			{
				{-1.f * surfaceScale, -1.f * surfaceScale, 1.f },
				{surfaceColor.x, surfaceColor.y, surfaceColor.z }
			},
			{
				{-1.f * surfaceScale, 1.f * surfaceScale, 1.f },
				{surfaceColor.x, surfaceColor.y, surfaceColor.z }
			}
	};

	std::vector<Vertex> pointsVertices;
	pointsVertices.reserve(gridSize.x * gridSize.y);
	glm::vec2 base = {-(gridSize.x - resolution) * resolution / 2.f, -(gridSize.y - resolution) * resolution / 2.f};
	for(unsigned int i = 0; i < gridSize.x; ++i)
	{
		for(unsigned j = 0; j < gridSize.y; ++j)
		{
			pointsVertices.push_back({
				{i * resolution + base.x, j * resolution + base.y, 1.f },
				{1.f, 1.f, 1.f }
			});
		}
	}
	/////////////////

	Renderable * surface = new Renderable(surfaceVertices, GL_TRIANGLES);
	surface->initVAO();
	Renderable * points = new Renderable{pointsVertices, GL_POINTS, pointSize};
	GLuint * vbo = points->getVBO();
	AcousticFDTD * fdtd = new AcousticFDTD(gridSize, vbo, blockSize, dataPerThread);
	points->initVAO();

	/////////////////

	auto appStart = std::chrono::high_resolution_clock::now();

	while(window.update())
	{
		auto appTime = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(appTime - appStart).count();

		shader.setUniform4m("projection", glm::value_ptr(projection));
		shader.setUniform4m("view", glm::value_ptr(view));
		shader.setUniform4m("model", glm::value_ptr(model));

		surface->draw();

		fdtd->draw();

		points->draw();
	}

	delete fdtd;
	delete points;
	delete surface;
	window.close();
	return 0;
}
