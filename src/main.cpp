#include <iostream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <ctime>

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
	unsigned int dataPerThread = 1;
	unsigned int blockSize = 32; // 32 x 32
	bool oglContext = false;
	for(unsigned int i = 1; i < argc; ++i)
	{
		if(argc > i + 1)
		{
			if (!strcmp(argv[i], "-r"))
			{
				resolution = 10.f/atof(argv[i + 1]);
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
		if (!strcmp(argv[i], "-ogl"))
		{
			oglContext = true;
		}
	}
	glm::ivec2 gridSize = {surfaceSize.x / resolution, surfaceSize.y / resolution};
	///////////////

	GLuint * vbo = nullptr;
	AcousticFDTD * fdtd = new AcousticFDTD(gridSize, vbo, blockSize, dataPerThread);
	while(true)
	{
		fdtd->draw();
	}

	return 0;
}
