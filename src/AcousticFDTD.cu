#include "AcousticFDTD.h"

#include <cmath>

const int CUDA_THREADS_X = 256;
const int CUDA_THREADS_Y = 256;

AcousticFDTD::AcousticFDTD(glm::ivec2 & gridSize, void * vertexPointer)
{
	_gridSize = gridSize;
	_vertexPointer = static_cast<glm::vec3 *>(vertexPointer);
	_cudaBlockSize = dim3(CUDA_THREADS_X, CUDA_THREADS_Y);
	const int bx = (gridSize.x + CUDA_THREADS_X - 1) / CUDA_THREADS_X;
	const int by = (gridSize.y + CUDA_THREADS_Y - 1) / CUDA_THREADS_Y;
	_cudaGridSize = dim3(bx, by);
	_dataPerThread = glm::vec2(1, 1);

	for(unsigned int i = 0; i < 2; ++i)
	{
		cudaMalloc((void **)&_grid[i], (gridSize.x + 1) * (gridSize.y + 1) * sizeof(SpacePoint));
		cudaMemset(_grid[i], 0, (gridSize.x + 1) * (gridSize.y + 1) * sizeof(SpacePoint));
		cudaMalloc((void **)&_murX[i], gridSize.y * 4 * sizeof(float));
		cudaMemset(_murX[i], 0, gridSize.y * 4 * sizeof(float));
		cudaMalloc((void **)&_murY[i], gridSize.x * 4 * sizeof(float));
		cudaMemset(_murY[i], 0, gridSize.x * 4 * sizeof(float));
	}
}

AcousticFDTD::~AcousticFDTD()
{
	cudaFree(_grid[0]); cudaFree(_grid[1]);
	cudaFree(_murX[0]); cudaFree(_murX[1]);
	cudaFree(_murY[0]); cudaFree(_murY[1]);
}

__global__ void updateV(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid,
		AcousticFDTD::SpacePoint * outGrid, float dtOverDx, float density)
{
	const int startI = blockIdx.x * blockDim.x + threadIdx.x;
	const int startJ = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int rangeI = (startI + dataPerThread.x) < gridSize.x ? (startI + dataPerThread.x) : gridSize.x;
	unsigned int rangeJ = (startJ + dataPerThread.y) < gridSize.y ? (startJ + dataPerThread.y) : gridSize.y;

	for(unsigned int i = startI == 0 ? 1 : startI; i < rangeI; ++i)
	{
		for(unsigned int j = startJ; j < rangeJ; ++j)
		{
			outGrid[((i * gridSize.y) + j)].vx += -dtOverDx / density * (inGrid[(i * gridSize.y) + j].soundPressure - inGrid[((i - 1) * gridSize.y) + j].soundPressure);
		}
	}

	for(unsigned int i = startI; i < rangeI; ++i)
	{
		for(unsigned int j = startJ == 0 ? 1 : startJ; j < rangeJ; ++j)
		{
			outGrid[((i * gridSize.y) + j)].vy += -dtOverDx / density * (inGrid[(i * gridSize.y) + j].soundPressure - inGrid[i * gridSize.y + (j - 1)].soundPressure);
		}
	}
}

__global__ void updateP(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid,
		AcousticFDTD::SpacePoint * outGrid, float dtOverDx, float bulkModulus)
{
	const int startI = blockIdx.x * blockDim.x + threadIdx.x;
	const int startJ = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int rangeI = (startI + dataPerThread.x) < gridSize.x ? (startI + dataPerThread.x) : gridSize.x;
	unsigned int rangeJ = (startJ + dataPerThread.y) < gridSize.y ? (startJ + dataPerThread.y) : gridSize.y;

	for(unsigned int i = startI; i < rangeI; ++i)
	{
		for(unsigned int j = startJ; j < rangeJ; ++j)
		{
			outGrid[((i * gridSize.y) + j)].soundPressure += -(bulkModulus * dtOverDx) * ((inGrid[((i + 1) * gridSize.y) + j].vx - inGrid[(i * gridSize.y) + j].vx) + (inGrid[(i * gridSize.y) + (j + 1)].vy - inGrid[(i * gridSize.y) + j].vy));
		}
	}
}

__global__ void mur2nd(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid, AcousticFDTD::SpacePoint * outGrid,
		float * murX[2], float * murY[2], float dt, float dx, float density, float bulkModulus)
{
	const int startI = blockIdx.x * blockDim.x + threadIdx.x;
	const int startJ = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int rangeI = (startI + dataPerThread.x) < (gridSize.x - 2) ? (startI + dataPerThread.x) : (gridSize.x - 2);
	unsigned int rangeJ = (startJ + dataPerThread.y) < (gridSize.y - 2) ? (startJ + dataPerThread.y) : (gridSize.y - 2);

	float v = sqrt(bulkModulus/density);	// Wave velocity
	int i,j;

	for(i = startI < 2 ? 2 : startI; i < rangeI; ++i){
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
	for(j = startJ < 2 ? 2 : startJ; j < rangeJ; ++j){
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

	// corners are computed by first thread:
	// Mur's 1st Order Absorption for 4 corners
	if(startI == 0)
	{
		i = 1;
		outGrid[i * gridSize.y].soundPressure = murY[0][i * 4 + 1] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + 1].soundPressure - murY[0][i * 4]);
		outGrid[i * gridSize.y + gridSize.y-1].soundPressure = murY[0][i * 4 + 2] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + gridSize.y-2].soundPressure - murY[0][i * 4 + 3]);
		i = gridSize.x-2;
		outGrid[i * gridSize.y].soundPressure = murY[0][i * 4 + 1] + (v*dt-dx)/(v*dt+dx) * (inGrid[i* gridSize.y + 1].soundPressure - murY[0][i * 4]);
		outGrid[i * gridSize.y + gridSize.y-1].soundPressure = murY[0][i * 4 + 2] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + gridSize.y-2].soundPressure - murY[0][i * 4 + 3]);
	}
	if(startJ == 0)
	{
		j = 1;
		outGrid[0 * gridSize.y + j].soundPressure = murX[0][1 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[1 * gridSize.y + j].soundPressure - murX[0][0 + j * 4]);
		outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = murX[0][2 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure - murX[0][3 + j * 4]);
		j = gridSize.y - 2;
		outGrid[0 * gridSize.y + j].soundPressure = murX[0][1 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[1 * gridSize.y + j].soundPressure - murX[0][0 + j * 4]);
		outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = murX[0][2+ j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure - murX[0][3 + j * 4]);
	}
}

__global__ void mur2ndCopy(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, float * murX[2], float * murY[2])
{
	const int startI = blockIdx.x * blockDim.x + threadIdx.x;
	const int startJ = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int rangeI = (startI + dataPerThread.x) < gridSize.x ? (startI + dataPerThread.x) : gridSize.x;
	unsigned int rangeJ = (startJ + dataPerThread.y) < gridSize.y ? (startJ + dataPerThread.y) : gridSize.y;

	for(unsigned int i = startI; i < rangeI; ++i)
	{
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
	for(unsigned int j = startJ; j < rangeJ; ++j){
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

__global__ void updateColors(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, glm::vec3 * vertexPointer)
{
	const int startI = blockIdx.x * blockDim.x + threadIdx.x;
	const int startJ = blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int rangeI = (startI + dataPerThread.x) < gridSize.x ? (startI + dataPerThread.x) : gridSize.x;
	unsigned int rangeJ = (startJ + dataPerThread.y) < gridSize.y ? (startJ + dataPerThread.y) : gridSize.y;

	//setGridColors
	for(unsigned int i = startI; i < rangeI; ++i)
	{
		for(unsigned int j = startJ; j < rangeJ; ++j)
		{
			float amplifier = 100.f;
			float grayScale = abs(grid[i * gridSize.y + j].soundPressure) * amplifier;
			// 2 * index + 1 because vbo consists of 2 vec3 (position and color)
			vertexPointer[2 * (i * gridSize.y + j) + 1] = {grayScale, grayScale, grayScale};
		}
	}
}

__global__ void updatePoint(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, glm::ivec2 point, float sigPoint)
{
	grid[point.x * gridSize.y + point.y].soundPressure += sigPoint;
}

void AcousticFDTD::draw()
{
	updateV<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)!_bufferSwap],
												_grid[(int)_bufferSwap], _dtOverDx, _density);
	updateP<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)!_bufferSwap],
												_grid[(int)_bufferSwap], _dtOverDx, _bulkModulus);

	mur2nd<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)!_bufferSwap],
												_grid[(int)_bufferSwap], _murX, _murY, _dt, _dx, _density, _bulkModulus);
	//copy previous values
	mur2ndCopy<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid, _murX, _murY);

	/* Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window) */
	if( _nPoint < (1.0/_freq)/_dt ){
		_sigPoint = (1.0-cos((2.0*M_PI*_freq*_nPoint*_dt)))/2.0 * sin((2.0*M_PI*_freq*_nPoint*_dt));
		updatePoint<<<1, 1>>>(_gridSize, _grid, _pointSource, _sigPoint);
	}

	updateColors<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid, _vertexPointer);

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