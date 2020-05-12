#include "AcousticFDTD.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <CudaHelper.h>

AcousticFDTD::AcousticFDTD(glm::ivec2 & gridSize, GLuint * vbo, unsigned int blockSize, unsigned int dataPerThread)
{
	_gridSize = gridSize;
	if(dataPerThread > 0 && dataPerThread <= gridSize.x)
	{
		_dataPerThread = glm::ivec2(dataPerThread, dataPerThread);
	}
	else
	{
		_dataPerThread = glm::ivec2(1, 1);
	}

	if(blockSize > 0 && blockSize <= 32)
	{
		_cudaBlockSize = dim3(blockSize, blockSize);
	}
	else
	{
		_cudaBlockSize = dim3(32, 32);
	}
	const int bx = (gridSize.x + _cudaBlockSize.x - 1) / _cudaBlockSize.x;
	const int by = (gridSize.y + _cudaBlockSize.y - 1) / _cudaBlockSize.y;
	_cudaGridSize = dim3(bx, by);

	CudaSafeCall(cudaGraphicsGLRegisterBuffer(&_cudaVboRes, *vbo, cudaGraphicsMapFlagsNone));
	
	for(unsigned int i = 0; i < 2; ++i)
	{
		CudaSafeCall(cudaMalloc((void **)&_grid[i], (gridSize.x + 1) * (gridSize.y + 1) * sizeof(SpacePoint)));
		CudaSafeCall(cudaMemset(_grid[i], 0, (gridSize.x + 1) * (gridSize.y + 1) * sizeof(SpacePoint)));
		CudaSafeCall(cudaMalloc((void **)&_murX[i], gridSize.y * 4 * sizeof(float)));
		CudaSafeCall(cudaMemset(_murX[i], 0, gridSize.y * 4 * sizeof(float)));
		CudaSafeCall(cudaMalloc((void **)&_murY[i], gridSize.x * 4 * sizeof(float)));
		CudaSafeCall(cudaMemset(_murY[i], 0, gridSize.x * 4 * sizeof(float)));
	}
}

AcousticFDTD::~AcousticFDTD()
{
	CudaSafeCall(cudaGraphicsUnregisterResource(_cudaVboRes));

	CudaSafeCall(cudaFree(_grid[0]));
	CudaSafeCall(cudaFree(_grid[1]));
	CudaSafeCall(cudaFree(_murX[0]));
	CudaSafeCall(cudaFree(_murX[1]));
	CudaSafeCall(cudaFree(_murY[0]));
	CudaSafeCall(cudaFree(_murY[1]));
}

__global__ void updateV(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * inGrid,
		AcousticFDTD::SpacePoint * outGrid, float dtOverDx, float density)
{
	const int startI = (blockIdx.x * blockDim.x + threadIdx.x) * dataPerThread.x;
	const int startJ = (blockIdx.y * blockDim.y + threadIdx.y) * dataPerThread.y;

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
	const int startI = (blockIdx.x * blockDim.x + threadIdx.x) * dataPerThread.x;
	const int startJ = (blockIdx.y * blockDim.y + threadIdx.y) * dataPerThread.y;

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
		float * murX1, float * murX2, float * murY1, float * murY2, float dt, float dx, float density, float bulkModulus)
{
	const int startI = (blockIdx.x * blockDim.x + threadIdx.x) * dataPerThread.x;
	const int startJ = (blockIdx.y * blockDim.y + threadIdx.y) * dataPerThread.y;

	unsigned int rangeI = (startI + dataPerThread.x) < (gridSize.x - 2) ? (startI + dataPerThread.x) : (gridSize.x - 2);
	unsigned int rangeJ = (startJ + dataPerThread.y) < (gridSize.y - 2) ? (startJ + dataPerThread.y) : (gridSize.y - 2);

	float v = sqrt(bulkModulus/density);	// Wave velocity
	int i,j;

	for(i = startI < 2 ? 2 : startI; i < rangeI; ++i){
		outGrid[(i * gridSize.y)].soundPressure = - murY2[(i * 4 + 1)]
		                                         + (v*dt-dx)/(v*dt+dx) * ( inGrid[(i * gridSize.y + 1)].soundPressure + murY2[(i * 4)] )
		                                         + (2.0*dx)/(v*dt+dx) * ( murY1[i * 4] + murY1[i * 4 + 1] )
		                                         + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                                           * ( murY1[(i+1) * 4] - 2.0 * murY1[i * 4]
		                                               + murY1[(i-1) * 4] + murY1[(i+1) * 4 + 1]
		                                               - 2.0 * murY1[i* 4 + 1] + murY1[(i-1) * 4 + 1] );
		outGrid[(i * gridSize.y) + (gridSize.y - 1)].soundPressure = - murY2[(i * 4) + 2]
		                                                             + (v*dt-dx)/(v*dt+dx) * ( inGrid[(i * gridSize.y) + gridSize.y-2].soundPressure + murY2[i * 4 + 3] )
		                                                             + (2.0*dx)/(v*dt+dx) * ( murY1[i * 4 + 3] + murY1[i * 4 + 2] )
		                                                             + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                                                               * ( murY1[(i+1) * 4 + 3] - 2.0 * murY1[i * 4 + 3]
		                                                                   + murY1[(i-1) * 4 + 3] + murY1[(i+1) * 4 + 2]
		                                                                   - 2.0 * murY1[i * 4 + 2] + murY1[(i-1) * 4 + 2] );
	}
	for(j = startJ < 2 ? 2 : startJ; j < rangeJ; ++j){
		outGrid[j].soundPressure = - murX2[1 + j * 4]
		                         + (v*dt-dx)/(v*dt+dx) * ( inGrid[1 * gridSize.y + j].soundPressure + murX2[0 + j * 4] )
		                         + (2.0*dx)/(v*dt+dx) * ( murX1[0 + j * 4] + murX1[1 + j * 4] )
		                         + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                           * ( murX1[0 + (j+1) * 4] - 2.0 * murX1[0 + j * 4]
		                               + murX1[0 + (j-1) * 4] + murX1[1 + (j+1) * 4]
		                               - 2.0 * murX1[1 + j * 4] + murX1[1 + (j-1) * 4] );
		outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = - murX2[2 + j * 4]
		                                                         + (v*dt-dx)/(v*dt+dx) * ( inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure + murX2[3 + j * 4] )
		                                                         + (2.0*dx)/(v*dt+dx) * ( murX1[3 + j * 4] + murX1[2 + j * 4] )
		                                                         + (dx*v*v*dt*dt)/(2.0*dx*dx*(v*dt+dx))
		                                                           * ( murX1[3 + (j+1) * 4] - 2.0 * murX1[3 + j * 4]
		                                                               + murX1[3 + (j-1) * 4] + murX1[2 + (j+1) * 4]
		                                                               - 2.0 * murX1[2 + j * 4] + murX1[2 + (j-1) * 4] );
	}

	// Mur's 1st Order Absorption for 4 corners
	// computed by first thread
	if(startI == 0)
	{
		i = 1;
		outGrid[i * gridSize.y].soundPressure = murY1[i * 4 + 1] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + 1].soundPressure - murY1[i * 4]);
		outGrid[i * gridSize.y + gridSize.y-1].soundPressure = murY1[i * 4 + 2] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + gridSize.y-2].soundPressure - murY1[i * 4 + 3]);
		i = gridSize.x-2;
		outGrid[i * gridSize.y].soundPressure = murY1[i * 4 + 1] + (v*dt-dx)/(v*dt+dx) * (inGrid[i* gridSize.y + 1].soundPressure - murY1[i * 4]);
		outGrid[i * gridSize.y + gridSize.y-1].soundPressure = murY1[i * 4 + 2] + (v*dt-dx)/(v*dt+dx) * (inGrid[i * gridSize.y + gridSize.y-2].soundPressure - murY1[i * 4 + 3]);
	}
	if(startJ == 0)
	{
		j = 1;
		outGrid[0 * gridSize.y + j].soundPressure = murX1[1 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[1 * gridSize.y + j].soundPressure - murX1[0 + j * 4]);
		outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = murX1[2 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure - murX1[3 + j * 4]);
		j = gridSize.y - 2;
		outGrid[0 * gridSize.y + j].soundPressure = murX1[1 + j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[1 * gridSize.y + j].soundPressure - murX1[0 + j * 4]);
		outGrid[(gridSize.x-1) * gridSize.y + j].soundPressure = murX1[2+ j * 4] + (v*dt-dx)/(v*dt+dx) * (inGrid[(gridSize.x-2) * gridSize.y + j].soundPressure - murX1[3 + j * 4]);
	}
}

__global__ void mur2ndCopy(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid,
		float * murX1, float * murX2, float * murY1, float * murY2)
{
	const int startI = (blockIdx.x * blockDim.x + threadIdx.x) * dataPerThread.x;
	const int startJ = (blockIdx.y * blockDim.y + threadIdx.y) * dataPerThread.y;

	unsigned int rangeI = (startI + dataPerThread.x) < gridSize.x ? (startI + dataPerThread.x) : gridSize.x;
	unsigned int rangeJ = (startJ + dataPerThread.y) < gridSize.y ? (startJ + dataPerThread.y) : gridSize.y;

	for(unsigned int i = startI; i < rangeI; ++i)
	{
		// Copy 1st Old Values to 2nd Old Values
		for(unsigned int j = 0; j < 4; ++j)
		{
			murY2[i * 4 + j] = murY1[i * 4 + j];
		}

		// Copy Present Values
		murY1[i * 4 + 0] = grid[i * gridSize.y + 0].soundPressure;
		murY1[i * 4 + 1] = grid[i * gridSize.y + 1].soundPressure;
		murY1[i * 4 + 2] = grid[i * gridSize.y + gridSize.y-2].soundPressure;
		murY1[i * 4 + 3] = grid[i * gridSize.y + gridSize.y-1].soundPressure;
	}
	for(unsigned int j = startJ; j < rangeJ; ++j){
		// Copy 1st Old Values to 2nd Old Values
		for(unsigned int i = 0; i < 4; ++i)
		{
			murX2[i + j * 4] = murX1[i + j * 4];
		}

		// Copy Present Values
		murX1[0 + j * 4] = grid[0 * gridSize.y + j].soundPressure;
		murX1[1 + j * 4] = grid[1 * gridSize.y + j].soundPressure;
		murX1[2 + j * 4] = grid[(gridSize.x-2) * gridSize.y + j].soundPressure;
		murX1[3 + j * 4] = grid[(gridSize.x-1) * gridSize.y + j].soundPressure;
	}
}

__global__ void updateColors(glm::ivec2 dataPerThread, glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, Vertex * vertexPointer)
{
	const int startI = (blockIdx.x * blockDim.x + threadIdx.x) * dataPerThread.x;
	const int startJ = (blockIdx.y * blockDim.y + threadIdx.y) * dataPerThread.y;

	unsigned int rangeI = (startI + dataPerThread.x) < gridSize.x ? (startI + dataPerThread.x) : gridSize.x;
	unsigned int rangeJ = (startJ + dataPerThread.y) < gridSize.y ? (startJ + dataPerThread.y) : gridSize.y;

	//setGridColors
	for(unsigned int i = startI; i < rangeI; ++i)
	{
		for(unsigned int j = startJ; j < rangeJ; ++j)
		{
			float amplifier = 100.f;
			float grayScale = abs(grid[i * gridSize.y + j].soundPressure) * amplifier;
			Vertex & v = vertexPointer[(i * gridSize.y + j)];
			v.color = { grayScale, 1.f, 1.f };
		}
	}
}

__global__ void updatePoint(glm::ivec2 gridSize, AcousticFDTD::SpacePoint * grid, glm::ivec2 point, float sigPoint)
{
	grid[point.x * gridSize.y + point.y].soundPressure += sigPoint;
}

void AcousticFDTD::draw()
{
	CudaSafeCall(cudaGraphicsMapResources(1, &_cudaVboRes, 0));
	size_t size;
	CudaSafeCall(cudaGraphicsResourceGetMappedPointer((void **)(&_vertexPointer), &size, _cudaVboRes));

	updateV<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)!_bufferSwap],
												_grid[(int)_bufferSwap], _dtOverDx, _density);
	cudaDeviceSynchronize();
	CudaCheckError();

	updateP<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)!_bufferSwap],
												_grid[(int)_bufferSwap], _dtOverDx, _bulkModulus);
	cudaDeviceSynchronize();
	CudaCheckError();

	mur2nd<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)_bufferSwap],
												_grid[(int)_bufferSwap], _murX[0], _murX[1], _murY[0], _murY[1],
												_dt, _dx, _density, _bulkModulus);
	cudaDeviceSynchronize();
	CudaCheckError();

	//copy previous values
	mur2ndCopy<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)_bufferSwap],
	                                            _murX[0], _murX[1], _murY[0], _murY[1]);
	cudaDeviceSynchronize();
	CudaCheckError();

	// Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window)
	if( _nPoint < (1.0/_freq)/_dt ){
		_sigPoint = (1.0-cos((2.0*M_PI*_freq*_nPoint*_dt)))/2.0 * sin((2.0*M_PI*_freq*_nPoint*_dt));
		updatePoint<<<1, 1>>>(_gridSize, _grid[(int)_bufferSwap], _pointSource, _sigPoint);
		cudaDeviceSynchronize();
		CudaCheckError();
	}

	updateColors<<<_cudaGridSize, _cudaBlockSize>>>(_dataPerThread, _gridSize, _grid[(int)_bufferSwap], _vertexPointer);
	cudaDeviceSynchronize();
	CudaCheckError();

	CudaSafeCall(cudaGraphicsUnmapResources(1, &_cudaVboRes, 0));

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
