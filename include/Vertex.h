#pragma once

#include <glad/glad.h>
#include <cuda.h>
#include <cuda_runtime.h>
#define GLM_FORCE_CUDA
#include "../dependencies/glm/glm/glm.hpp"

struct Vertex
{
	glm::vec4 position;
	glm::vec4 color;
};