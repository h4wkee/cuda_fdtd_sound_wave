#include "AcousticFDTD.h"

#include <omp.h>

#define _USE_MATH_DEFINES
#define FDTD_LOOP _Pragma("omp for collapse(2)")

#include <chrono>

AcousticFDTD::AcousticFDTD(glm::vec3 * gridColors, glm::ivec2 & gridSize)
{
	omp_set_nested(1);
	_gridColors = gridColors;

	_gridSize = gridSize;
	_grid = new SpacePoint[(gridSize.x + 1) * (gridSize.y + 1)]{};
	_murX1 = new float[gridSize.y * 4]{};
	_murX2 = new float[gridSize.y * 4]{};
	_murY1 = new float[gridSize.x * 4]{};
	_murY2 = new float[gridSize.x * 4]{};
}

AcousticFDTD::~AcousticFDTD()
{
	delete[] _grid;
	delete[] _murX1;
	delete[] _murX2;
	delete[] _murY1;
	delete[] _murY2;
}

void AcousticFDTD::updateV()
{
	FDTD_LOOP
	for(unsigned int i = 1; i < _gridSize.x; ++i)
	{
		for(unsigned int j = 0; j < _gridSize.y; ++j)
		{
//			int threadid=omp_get_thread_num();
//			printf("Thread: %d: i=%d j=%d\n", threadid, i, j);
			_grid[((i * _gridSize.y) + j)].vx += -_dtOverDx / _density * (_grid[(i * _gridSize.y) + j].soundPressure - _grid[((i - 1) * _gridSize.y) + j].soundPressure);
		}
	}

	FDTD_LOOP
	for(unsigned int i = 0; i < _gridSize.x; ++i)
	{
		for(unsigned int j = 1; j < _gridSize.y; ++j)
		{
			_grid[((i * _gridSize.y) + j)].vy += -_dtOverDx / _density * (_grid[(i * _gridSize.y) + j].soundPressure - _grid[i * _gridSize.y + (j - 1)].soundPressure);
		}
	}
}

void AcousticFDTD::updateP()
{
	FDTD_LOOP
	for(unsigned int i = 0; i < _gridSize.x; ++i)
	{
		for(unsigned int j = 0; j < _gridSize.y; ++j)
		{
			_grid[((i * _gridSize.y) + j)].soundPressure += -(_bulkModulus * _dtOverDx) * ((_grid[((i + 1) * _gridSize.y) + j].vx - _grid[(i * _gridSize.y) + j].vx) + (_grid[(i * _gridSize.y) + (j + 1)].vy - _grid[(i * _gridSize.y) + j].vy));
		}
	}
}

/* Mur's 2nd Order Absorption */
void AcousticFDTD::mur2nd()
{
	#pragma omp for
	for(int i=2;i<_gridSize.x-2;i++){
		_grid[(i * _gridSize.y)].soundPressure = - _murY2[(i * 4 + 1)]
			+ (_v*_dt-_dx)/(_v*_dt+_dx) * ( _grid[(i * _gridSize.y + 1)].soundPressure + _murY2[(i * 4)] )
            + (2.0*_dx)/(_v*_dt+_dx) * ( _murY1[i * 4] + _murY1[i * 4 + 1] )
            + (_dx*_v*_v*_dt*_dt)/(2.0*_dx*_dx*(_v*_dt+_dx))
            * ( _murY1[(i+1) * 4] - 2.0 * _murY1[i * 4]
                + _murY1[(i-1) * 4] + _murY1[(i+1) * 4 + 1]
                - 2.0 * _murY1[i* 4 + 1] + _murY1[(i-1) * 4 + 1] );
		_grid[(i * _gridSize.y) + (_gridSize.y - 1)].soundPressure = - _murY2[(i * 4) + 2]
             + (_v*_dt-_dx)/(_v*_dt+_dx) * ( _grid[(i * _gridSize.y) + _gridSize.y-2].soundPressure + _murY2[i * 4 + 3] )
             + (2.0*_dx)/(_v*_dt+_dx) * ( _murY1[i * 4 + 3] + _murY1[i * 4 + 2] )
             + (_dx*_v*_v*_dt*_dt)/(2.0*_dx*_dx*(_v*_dt+_dx))
             * ( _murY1[(i+1) * 4 + 3] - 2.0 * _murY1[i * 4 + 3]
                + _murY1[(i-1) * 4 + 3] + _murY1[(i+1) * 4 + 2]
                - 2.0 * _murY1[i * 4 + 2] + _murY1[(i-1) * 4 + 2] );
	}
	#pragma omp for
	for(int j=2;j<_gridSize.y-2;j++){
		_grid[j].soundPressure = - _murX2[1 + j * 4]
			+ (_v*_dt-_dx)/(_v*_dt+_dx) * ( _grid[1 * _gridSize.y + j].soundPressure + _murX2[0 + j * 4] )
			+ (2.0*_dx)/(_v*_dt+_dx) * ( _murX1[0 + j * 4] + _murX1[1 + j * 4] )
			+ (_dx*_v*_v*_dt*_dt)/(2.0*_dx*_dx*(_v*_dt+_dx))
			* ( _murX1[0 + (j+1) * 4] - 2.0 * _murX1[0 + j * 4]
				+ _murX1[0 + (j-1) * 4] + _murX1[1 + (j+1) * 4]
				- 2.0 * _murX1[1 + j * 4] + _murX1[1 + (j-1) * 4] );
		_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = - _murX2[2 + j * 4]
	         + (_v*_dt-_dx)/(_v*_dt+_dx) * ( _grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure + _murX2[3 + j * 4] )
	         + (2.0*_dx)/(_v*_dt+_dx) * ( _murX1[3 + j * 4] + _murX1[2 + j * 4] )
	         + (_dx*_v*_v*_dt*_dt)/(2.0*_dx*_dx*(_v*_dt+_dx))
	           * ( _murX1[3 + (j+1) * 4] - 2.0 * _murX1[3 + j * 4]
	                 + _murX1[3 + (j-1) * 4] + _murX1[2 + (j+1) * 4]
	                 - 2.0 * _murX1[2 + j * 4] + _murX1[2 + (j-1) * 4] );
	}

	#pragma omp single
	{
		/* Mur's 1st Order Absorption for 4 corners*/
		int i = 1;
		_grid[i * _gridSize.y].soundPressure = _murY1[i * 4 + 1] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[i * _gridSize.y + 1].soundPressure - _murY1[i * 4]);
		_grid[i * _gridSize.y + _gridSize.y-1].soundPressure = _murY1[i * 4 + 2] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[i * _gridSize.y + _gridSize.y-2].soundPressure - _murY1[i * 4 + 3]);
		i = _gridSize.x-2;
		_grid[i * _gridSize.y].soundPressure = _murY1[i * 4 + 1] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[i* _gridSize.y + 1].soundPressure - _murY1[i * 4]);
		_grid[i * _gridSize.y + _gridSize.y-1].soundPressure = _murY1[i * 4 + 2] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[i * _gridSize.y + _gridSize.y-2].soundPressure - _murY1[i * 4 + 3]);
		int j = 1;
		_grid[0 * _gridSize.y + j].soundPressure = _murX1[1 + j * 4] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[1 * _gridSize.y + j].soundPressure - _murX1[0 + j * 4]);
		_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = _murX1[2 + j * 4] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure - _murX1[3 + j * 4]);
		j = _gridSize.y - 2;
		_grid[0 * _gridSize.y + j].soundPressure = _murX1[1 + j * 4] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[1 * _gridSize.y + j].soundPressure - _murX1[0 + j * 4]);
		_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = _murX1[2+ j * 4] + (_v*_dt-_dx)/(_v*_dt+_dx) * (_grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure - _murX1[3 + j * 4]);
	};
}

/* Copy the Filed Values for Mur's 2nd Order Absorption */
void AcousticFDTD::mur2ndCopy()
{
	#pragma omp for
	for(int i=0;i<_gridSize.x;i++){
		/* Copy 1st Old Values to 2nd Old Values*/
		_murY2[i * 4 + 0] = _murY1[i * 4 + 0];
		_murY2[i * 4 + 1] = _murY1[i * 4 + 1];
		_murY2[i * 4 + 2] = _murY1[i * 4 + 2];
		_murY2[i * 4 + 3] = _murY1[i * 4 + 3];

		/* Copy Present Values */
		_murY1[i * 4 + 0] = _grid[i * _gridSize.y + 0].soundPressure;
		_murY1[i * 4 + 1] = _grid[i * _gridSize.y + 1].soundPressure;
		_murY1[i * 4 + 2] = _grid[i * _gridSize.y + _gridSize.y-2].soundPressure;
		_murY1[i * 4 + 3] = _grid[i * _gridSize.y + _gridSize.y-1].soundPressure;
	}
	#pragma omp for
	for(int j=0;j<_gridSize.y;j++){
		/* Copy 1st Old Values to 2nd Old Values*/
		_murX2[0 + j * 4] = _murX1[0 + j * 4];
		_murX2[1 + j * 4] = _murX1[1 + j * 4];
		_murX2[2 + j * 4] = _murX1[2 + j * 4];
		_murX2[3 + j * 4] = _murX1[3 + j * 4];

		/* Copy Present Values */
		_murX1[0 + j * 4] = _grid[0 * _gridSize.y + j].soundPressure;
		_murX1[1 + j * 4] = _grid[1 * _gridSize.y + j].soundPressure;
		_murX1[2 + j * 4] = _grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure;
		_murX1[3 + j * 4] = _grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure;
	}
}

void AcousticFDTD::draw()
{
	auto loopStart = std::chrono::high_resolution_clock::now();

	#pragma omp parallel default(none)
	{
		updateV();

		# pragma omp barrier

		updateP();

		# pragma omp barrier

		mur2nd();

		# pragma omp barrier

		mur2ndCopy();

		# pragma omp single
		{
			/* Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window) */
			if (_nPoint < (1.0 / _freq) / _dt)
			{
				_sigPoint = (1.0 - cos((2.0 * M_PI * _freq * _nPoint * _dt))) / 2.0 *
				            sin((2.0 * M_PI * _freq * _nPoint * _dt));
				_grid[_pointSource.x * _gridSize.y + _pointSource.y].soundPressure += _sigPoint;
			}
		}
		//setGridColors
		FDTD_LOOP
		for(unsigned int i = 0; i < _gridSize.x; ++i)
		{
			for(unsigned int j = 0; j < _gridSize.y; ++j)
			{
				float amplifier = 100.f;
				float grayScale = abs(_grid[i * _gridSize.y + j].soundPressure) * amplifier;
				//std::cout << "GRAY SCALE: " << grayScale << std::endl;
				_gridColors[i * _gridSize.y + j] = {grayScale, grayScale, grayScale};
			}
		}

		#pragma omp barrier
	}
	auto loopEnd = std::chrono::high_resolution_clock::now();

	if(_nPoint % _randomPointSourceInterval == 0)
	{
		auto loopTime = std::chrono::duration_cast<std::chrono::microseconds>(loopEnd - loopStart).count();
		std::cout << "Loop time: " << loopTime / 1e3 << "ms\n";
		_nPoint = 0;
		randomPointSource();
		randomPointSource();
	}
	++_nPoint;

}

void AcousticFDTD::randomPointSource()
{
	_pointSource.x = rand() % _gridSize.x;
	_pointSource.y = rand() % _gridSize.y;
	_nPoint = 0;
}
