#pragma once

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

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
	AcousticFDTD(glm::vec3 * gridColors, glm::ivec2 & gridSize);
	~AcousticFDTD();

	void updateV();
	void updateP();
	void mur2nd();
	void mur2ndCopy();

	void draw();
	void randomPointSource();
private:

	glm::ivec2 _gridSize = {};
	glm::ivec2 _pointSource = {};
	int _randomPointSourceInterval = 250;
	int _nPoint = 0;
	float _sigPoint = 0.f;

	glm::vec3 * _gridColors;

	SpacePoint * _grid;
	float * _murX1;
	float * _murX2;
	float * _murY1;
	float * _murY2;

	float _dx = 10.0e-3;    // Spatial Resolution [m/space_point]
	float _dt = 15.0e-6;    // Temporal Resolution [s/step]
	float _dtOverDx = _dt/_dx;

	float _density   = 1.29;    // [kg/m^3]
	float _bulkModulus = 142.0e3;	// [Pa]
	float _freq = 1.0e3;    // Frequency of Initial Waveform [Hz]
	float _v = sqrt(_bulkModulus/_density);	// Wave velocity
};
