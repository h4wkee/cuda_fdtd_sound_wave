#define _USE_MATH_DEFINES
#include <iostream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <chrono>

#include <AcousticFDTD.h>

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
	}
	glm::ivec2 gridSize = {surfaceSize.x / resolution, surfaceSize.y / resolution};
	glm::vec3 * gridColors = new glm::vec3[gridSize.x * gridSize.y];
	AcousticFDTD fdtd(nullptr, gridSize);
	///////////////

	while(true)
	{
		fdtd.draw();
	}

	return 0;
}
