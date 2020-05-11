#pragma once

#include <iostream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
//#include <helper_cuda.h>
//#include <helper_cuda_gl.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/vector_angle.hpp>

class AcousticFDTD
{
public:
	struct SpacePoint
	{
		float vx;
		float vy;
		float soundPressure;
	};
	AcousticFDTD(glm::ivec2 & gridSize, GLuint * vbo);
	~AcousticFDTD();

	void draw();
	void randomPointSource();
private:

	glm::ivec2 _gridSize = {};
	glm::ivec2 _pointSource = {};
	int _randomPointSourceInterval = 100;
	int _nPoint = 0;
	float _sigPoint = 0.f;

	SpacePoint * _grid[2];
	float * _murX[2];
	float * _murY[2];

	float _dx = 10.0e-3;    // Spatial Resolution [m/space_point]
	float _dt = 15.0e-6;    // Temporal Resolution [s/step]
	float _dtOverDx = _dt/_dx;

	float _density   = 1.29;    // [kg/m^3]
	float _bulkModulus = 142.0e3;	// [Pa]
	float _freq = 1.0e3;    // Frequency of Initial Waveform [Hz]

	glm::vec3 * _vertexPointer;
	struct cudaGraphicsResource * _cudaVboRes;
	//GLuint _vbo;
	bool _bufferSwap = 0;

	dim3 _cudaBlockSize;
	dim3 _cudaGridSize;
	glm::ivec2 _dataPerThread;
};
