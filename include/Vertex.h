#pragma once

#include <glad/glad.h>
#include <cuda.h>
#include <cuda_runtime.h>
#define GLM_FORCE_CUDA
//#include "../dependencies/glm/glm/glm.hpp"
#include <glm/glm.hpp>

struct Vertex
{
	glm::vec3 position;
	glm::vec3 color;
	int padding[2]; // to 32 bytes
	//int padding[10]; // to 64 bytes
};