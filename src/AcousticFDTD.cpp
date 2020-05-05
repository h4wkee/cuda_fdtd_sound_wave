#include "AcousticFDTD.h"

#import <cmath>

AcousticFDTD::AcousticFDTD(glm::ivec2 & gridSize, void * vertexPointer)
{
	_gridSize = gridSize;
	_vertexPointer = vertexPointer;
//	_grid = new SpacePoint[(gridSize.x + 1) * (gridSize.y + 1)]{};
//	_murX[0] = new float[gridSize.y * 4]{};
//	_murX[1] = new float[gridSize.y * 4]{};
//	_murY[0] = new float[gridSize.x * 4]{};
//	_murY[1] = new float[gridSize.x * 4]{};

	for(unsigned int i = 0; i < 2; ++i)
	{
		cudaMalloc((void **)&_grid[i], (gridSize.x + 1) * (gridSize.y + 1) * sizeof(SpacePoint));
		cudaMemset(_grid[i], 0, (gridSize.x + 1) * (gridSize.y + 1) * sizeof(SpacePoint));
		cudaMalloc((void **)&_murX[0][i], gridSize.y * 4 * sizeof(float));
		cudaMemset(_murX[0][i], 0, gridSize.y * 4 * sizeof(float));
		cudaMalloc((void **)&_murX[1][i], gridSize.y * 4 * sizeof(float));
		cudaMemset(_murX[1][i], 0, gridSize.y * 4 * sizeof(float));
		cudaMalloc((void **)&_murY[0][i], gridSize.x * 4 * sizeof(float));
		cudaMemset(_murY[0][i], 0, gridSize.x * 4 * sizeof(float));
		cudaMalloc((void **)&_murY[1][i], gridSize.x * 4 * sizeof(float));
		cudaMemset(_murY[1][i], 0, gridSize.x * 4 * sizeof(float));
	}
}

AcousticFDTD::~AcousticFDTD()
{
//	delete[] _grid;
//	delete[] _murX[0];
//	delete[] _murX[1];
//	delete[] _murY[0];
//	delete[] _murY[1];

	cudaFree(_grid[0]); cudaFree(_grid[1]);
	cudaFree(_murX[0][0]); cudaFree(_murX[0][1]);
	cudaFree(_murX[1][0]); cudaFree(_murX[1][1]);
	cudaFree(_murY[0][0]); cudaFree(_murY[0][1]);
	cudaFree(_murY[1][0]); cudaFree(_murY[1][1]);
}

__global__ void updateV(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid, AcousticFDTD::SpacePoint * outGrid, float dtOverDx, float density)
{
	for(unsigned int i = 1; i < gridSize.x; ++i)
	{
		for(unsigned int j = 0; j < gridSize.y; ++j)
		{
			outGrid[((i * gridSize.y) + j)].vx += -dtOverDx / density * (inGrid[(i * gridSize.y) + j].soundPressure - inGrid[((i - 1) * gridSize.y) + j].soundPressure);
		}
	}

	for(unsigned int i = 0; i < gridSize.x; ++i)
	{
		for(unsigned int j = 1; j < gridSize.y; ++j)
		{
			outGrid[((i * gridSize.y) + j)].vy += -dtOverDx / density * (inGrid[(i * gridSize.y) + j].soundPressure - inGrid[i * gridSize.y + (j - 1)].soundPressure);
		}
	}
}

__global__ void updateP(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid, AcousticFDTD::SpacePoint * outGrid, float dtOverDx, float bulkModulus)
{
	for(unsigned int i = 0; i < gridSize.x; ++i)
	{
		for(unsigned int j = 0; j < gridSize.y; ++j)
		{
			outGrid[((i * gridSize.y) + j)].soundPressure += -(bulkModulus * dtOverDx) * ((inGrid[((i + 1) * gridSize.y) + j].vx - inGrid[(i * gridSize.y) + j].vx) + (inGrid[(i * gridSize.y) + (j + 1)].vy - inGrid[(i * gridSize.y) + j].vy));
		}
	}
}

__global__ void mur2nd(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid, AcousticFDTD::SpacePoint * outGrid,
		float * murX[2], float * murY[2], float dt, float dx, float density, float bulkModulus)
{
	float v = sqrt(bulkModulus/density);	// Wave velocity
	int i,j;

	for(i=2;i<gridSize.x-2;i++){
		outGrid[(i * gridSize.y)].soundPressure = - murY[1][(i * 4 + 1)]
		                                         + (v*dt-dx)/(v*dt+dx) * ( inGrid[(i * gridSize.y + 1)].soundPressure + murY[1][(i * 4)] )
		                                         + (2.0*dx)/(v*dt+dx) * ( murY[0][i * 4] + murY[0][i * 4 + 1] )
		                                         + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                                           * ( murY[0][(i+1) * 4] - 2.0 * murY[0][i * 4]
		                                               + murY[0][(i-1) * 4] + murY[0][(i+1) * 4 + 1]
		                                               - 2.0 * murY[0][i* 4 + 1] + murY[0][(i-1) * 4 + 1] );
		outGrid[(i * gridSize.y) + (gridSize.y - 1)].soundPressure = - murY[1][(i * 4) + 2]
		                                                             + (v*dt-dx)/(v*dt+dx) * ( inGrid[(i * gridSize.y) + gridSize.y-2].soundPressure + murY[1][i * 4 + 3] )
		                                                             + (2.0*dx)/(v*dt+dx) * ( murY[0][i * 4 + 3] + murY[0][i * 4 + 2] )
		                                                             + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                                                               * ( murY[0][(i+1) * 4 + 3] - 2.0 * murY[0][i * 4 + 3]
		                                                                   + murY[0][(i-1) * 4 + 3] + murY[0][(i+1) * 4 + 2]
		                                                                   - 2.0 * murY[0][i * 4 + 2] + murY[0][(i-1) * 4 + 2] );
	}
	for(j=2;j<gridSize.y-2;j++){
		outGrid[j].soundPressure = - murX[1][1 + j * 4]
		                         + (v*dt-dx)/(v*dt+dx) * ( inGrid[1 * gridSize.y + j].soundPressure + murX[1][0 + j * 4] )
		                         + (2.0*dx)/(v*dt+dx) * ( murX[0][0 + j * 4] + murX[0][1 + j * 4] )
		                         + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                           * ( murX[0][0 + (j+1) * 4] - 2.0 * murX[0][0 + j * 4]
		                               + murX[0][0 + (j-1) * 4] + murX[0][1 + (j+1) * 4]
		                               - 2.0 * murX[0][1 + j * 4] + murX[0][1 + (j-1) * 4] );
		outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = - murX[1][2 + j * 4]
		                                                         + (v*dt-dx)/(v*dt+dx) * ( inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure + murX[1][3 + j * 4] )
		                                                         + (2.0*dx)/(v*dt+dx) * ( murX[0][3 + j * 4] + murX[0][2 + j * 4] )
		                                                         + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                                                           * ( murX[0][3 + (j+1) * 4] - 2.0 * murX[0][3 + j * 4]
		                                                               + murX[0][3 + (j-1) * 4] + murX[0][2 + (j+1) * 4]
		                                                               - 2.0 * murX[0][2 + j * 4] + murX[0][2 + (j-1) * 4] );
	}

	/* Mur's 1st Order Absorption for 4 corners*/
	i = 1;
	outGrid[i * gridSize.y].soundPressure = murY[0][i * 4 + 1] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + 1].soundPressure - murY[0][i * 4]);
	outGrid[i * gridSize.y + gridSize.y-1].soundPressure = murY[0][i * 4 + 2] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + gridSize.y-2].soundPressure - murY[0][i * 4 + 3]);
	i = gridSize.x-2;
	outGrid[i * gridSize.y].soundPressure = murY[0][i * 4 + 1] + (v*dt-dx)/(v*dt+dx) * (inGrid[i* gridSize.y + 1].soundPressure - murY[0][i * 4]);
	outGrid[i * gridSize.y + gridSize.y-1].soundPressure = murY[0][i * 4 + 2] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + gridSize.y-2].soundPressure - murY[0][i * 4 + 3]);
	j = 1;
	outGrid[0 * gridSize.y + j].soundPressure = murX[0][1 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[1 * gridSize.y + j].soundPressure - murX[0][0 + j * 4]);
	outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = murX[0][2 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure - murX[0][3 + j * 4]);
	j = gridSize.y - 2;
	outGrid[0 * gridSize.y + j].soundPressure = murX[0][1 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[1 * gridSize.y + j].soundPressure - murX[0][0 + j * 4]);
	outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = murX[0][2+ j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure - murX[0][3 + j * 4]);
}

__global__ void mur2ndCopy(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, float * murX[2], float * murY[2])
{
	for(int i=0;i<gridSize.x;i++){
		/* Copy 1st Old Values to 2nd Old Values*/
		for(unsigned int i = 0; i < 4; ++i)
		{
			murY[1][i * 4 + i] = murY[0][i * 4 + i];
		}

		/* Copy Present Values */
		murY[0][i * 4 + 0] = grid[i * gridSize.y + 0].soundPressure;
		murY[0][i * 4 + 1] = grid[i * gridSize.y + 1].soundPressure;
		murY[0][i * 4 + 2] = grid[i * gridSize.y + gridSize.y-2].soundPressure;
		murY[0][i * 4 + 3] = grid[i * gridSize.y + gridSize.y-1].soundPressure;
	}
	for(int j=0;j<gridSize.y;j++){
		/* Copy 1st Old Values to 2nd Old Values*/
		for(unsigned int i = 0; i < 4; ++i)
		{
			murX[1][i + j * 4] = murX[0][i + j * 4];
		}

		/* Copy Present Values */
		murX[0][0 + j * 4] = grid[0 * gridSize.y + j].soundPressure;
		murX[0][1 + j * 4] = grid[1 * gridSize.y + j].soundPressure;
		murX[0][2 + j * 4] = grid[(gridSize.x-2) * gridSize.y + j].soundPressure;
		murX[0][3 + j * 4] = grid[(gridSize.x-1) * gridSize.y + j].soundPressure;
	}
}

__global__ void updateColors(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, void * vertexPointer)
{
	//setGridColors
	for(unsigned int i = 0; i < gridSize.x; ++i)
	{
		for(unsigned int j = 0; j < gridSize.y; ++j)
		{
			float amplifier = 100.f;
			float grayScale = abs(grid[i * gridSize.y + j].soundPressure) * amplifier;
			//std::cout << "GRAY SCALE: " << grayScale << std::endl;
			//_gridColors[i * _gridSize.y + j] = {grayScale, grayScale, grayScale};
			//TODO: assign to vertexPointer appropriate color
		}
	}
}
void AcousticFDTD::draw()
{
	updateV();
	updateP();

	mur2nd();
	//copy previous values
	mur2ndCopy();

	/* Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window) */
	if( _nPoint < (1.0/_freq)/_dt ){
		_sigPoint = (1.0-cos((2.0*M_PI*_freq*_nPoint*_dt)))/2.0 * sin((2.0*M_PI*_freq*_nPoint*_dt));
		//TODO: cudaMemCpy to cpu?? or run on gpu??
		_grid[_pointSource.x * _gridSize.y + _pointSource.y].soundPressure += _sigPoint;
	}

	updateColors(_gridSize, _grid, _vertexPointer);

	++_nPoint;

	if(_nPoint > _randomPointSourceInterval)
	{
		randomPointSource();
	}

	_bufferSwap = !_bufferSwap;
}

void AcousticFDTD::randomPointSource()
{
	_pointSource.x = rand() % _gridSize.x;
	_pointSource.y = rand() % _gridSize.y;
	_nPoint = 0;
}
