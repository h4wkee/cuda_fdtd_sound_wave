#include "AcousticFDTD.h"

#define _USE_MATH_DEFINES
#include <cmath>

AcousticFDTD::AcousticFDTD(glm::vec3 * gridColors, glm::ivec2 & gridSize)
{
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
	for(unsigned int i = 1; i < _gridSize.x; ++i)
	{
		for(unsigned int j = 0; j < _gridSize.y; ++j)
		{
			_grid[((i * _gridSize.y) + j)].vx += -_dtOverDx / _density * (_grid[(i * _gridSize.y) + j].soundPressure - _grid[((i - 1) * _gridSize.y) + j].soundPressure);
		}
	}
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
	float v = sqrt(_bulkModulus/_density);	// Wave velocity
	int i,j;

	for(i=2;i<_gridSize.x-2;i++){
		_grid[(i * _gridSize.y)].soundPressure = - _murY2[(i * 4 + 1)]
			+ (v*_dt-_dx)/(v*_dt+_dx) * ( _grid[(i * _gridSize.y + 1)].soundPressure + _murY2[(i * 4)] )
            + (2.0*_dx)/(v*_dt+_dx) * ( _murY1[i * 4] + _murY1[i * 4 + 1] )
            + (_dx*v*v*_dt*_dt)/(2.0*_dx*_dx*(v*_dt+_dx))
            * ( _murY1[(i+1) * 4] - 2.0 * _murY1[i * 4]
                + _murY1[(i-1) * 4] + _murY1[(i+1) * 4 + 1]
                - 2.0 * _murY1[i* 4 + 1] + _murY1[(i-1) * 4 + 1] );
		_grid[(i * _gridSize.y) + (_gridSize.y - 1)].soundPressure = - _murY2[(i * 4) + 2]
             + (v*_dt-_dx)/(v*_dt+_dx) * ( _grid[(i * _gridSize.y) + _gridSize.y-2].soundPressure + _murY2[i * 4 + 3] )
             + (2.0*_dx)/(v*_dt+_dx) * ( _murY1[i * 4 + 3] + _murY1[i * 4 + 2] )
             + (_dx*v*v*_dt*_dt)/(2.0*_dx*_dx*(v*_dt+_dx))
             * ( _murY1[(i+1) * 4 + 3] - 2.0 * _murY1[i * 4 + 3]
                + _murY1[(i-1) * 4 + 3] + _murY1[(i+1) * 4 + 2]
                - 2.0 * _murY1[i * 4 + 2] + _murY1[(i-1) * 4 + 2] );
	}
	for(j=2;j<_gridSize.y-2;j++){
		_grid[j].soundPressure = - _murX2[1 + j * 4]
			+ (v*_dt-_dx)/(v*_dt+_dx) * ( _grid[1 * _gridSize.y + j].soundPressure + _murX2[0 + j * 4] )
			+ (2.0*_dx)/(v*_dt+_dx) * ( _murX1[0 + j * 4] + _murX1[1 + j * 4] )
			+ (_dx*v*v*_dt*_dt)/(2.0*_dx*_dx*(v*_dt+_dx))
			* ( _murX1[0 + (j+1) * 4] - 2.0 * _murX1[0 + j * 4]
				+ _murX1[0 + (j-1) * 4] + _murX1[1 + (j+1) * 4]
				- 2.0 * _murX1[1 + j * 4] + _murX1[1 + (j-1) * 4] );
		_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = - _murX2[2 + j * 4]
	         + (v*_dt-_dx)/(v*_dt+_dx) * ( _grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure + _murX2[3 + j * 4] )
	         + (2.0*_dx)/(v*_dt+_dx) * ( _murX1[3 + j * 4] + _murX1[2 + j * 4] )
	         + (_dx*v*v*_dt*_dt)/(2.0*_dx*_dx*(v*_dt+_dx))
	           * ( _murX1[3 + (j+1) * 4] - 2.0 * _murX1[3 + j * 4]
	                 + _murX1[3 + (j-1) * 4] + _murX1[2 + (j+1) * 4]
	                 - 2.0 * _murX1[2 + j * 4] + _murX1[2 + (j-1) * 4] );
	}

	/* Mur's 1st Order Absorption for 4 corners*/
	i = 1;
	_grid[i * _gridSize.y].soundPressure = _murY1[i * 4 + 1] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[i * _gridSize.y + 1].soundPressure - _murY1[i * 4]);
	_grid[i * _gridSize.y + _gridSize.y-1].soundPressure = _murY1[i * 4 + 2] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[i * _gridSize.y + _gridSize.y-2].soundPressure - _murY1[i * 4 + 3]);
	i = _gridSize.x-2;
	_grid[i * _gridSize.y].soundPressure = _murY1[i * 4 + 1] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[i* _gridSize.y + 1].soundPressure - _murY1[i * 4]);
	_grid[i * _gridSize.y + _gridSize.y-1].soundPressure = _murY1[i * 4 + 2] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[i * _gridSize.y + _gridSize.y-2].soundPressure - _murY1[i * 4 + 3]);
	j = 1;
	_grid[0 * _gridSize.y + j].soundPressure = _murX1[1 + j * 4] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[1 * _gridSize.y + j].soundPressure - _murX1[0 + j * 4]);
	_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = _murX1[2 + j * 4] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure - _murX1[3 + j * 4]);
	j = _gridSize.y - 2;
	_grid[0 * _gridSize.y + j].soundPressure = _murX1[1 + j * 4] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[1 * _gridSize.y + j].soundPressure - _murX1[0 + j * 4]);
	_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = _murX1[2+ j * 4] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure - _murX1[3 + j * 4]);

	/* Copy Previous Values */
	mur2ndCopy();
}

/* Copy the Filed Values for Mur's 2nd Order Absorption */
void AcousticFDTD::mur2ndCopy()
{

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

void AcousticFDTD::mur1st()
{
	float v = sqrt(_bulkModulus/_density);

	for(int i=1;i<_gridSize.x-1;i++){
		_grid[i * _gridSize.y + 0].soundPressure = _murY1[i * 4 + 1] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[i * _gridSize.y + 1].soundPressure - _murY1[i * 4 + 0]);
		_grid[i * _gridSize.y + _gridSize.y-1].soundPressure = _murY1[i * 4 + 2] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[i * _gridSize.y + _gridSize.y-2].soundPressure - _murY1[i * 4 + 3]);
	}
	for(int j=1;j<_gridSize.y-1;j++){
		_grid[0 * _gridSize.y + j].soundPressure = _murX1[1 + j * 4] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[1 * _gridSize.y + j].soundPressure - _murX1[0 + j * 4]);
		_grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure = _murX1[2 + j * 4] + (v*_dt-_dx)/(v*_dt+_dx) * (_grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure - _murX1[3 + j * 4]);
	}

	/* Copy Previous Values */
	mur1stCopy();
}

/* Copy the Filed Values for Mur's 1st Order Absorption */
void AcousticFDTD::mur1stCopy()
{
	/* Copy Previous Values */
	for(int i=0;i<_gridSize.x;i++){
		_murY1[i * 4 + 0] = _grid[i * _gridSize.y + 0].soundPressure;
		_murY1[i * 4 + 1] = _grid[i * _gridSize.y + 1].soundPressure;
		_murY1[i * 4 + 2] = _grid[i * _gridSize.y + _gridSize.y-2].soundPressure;
		_murY1[i * 4 + 3] = _grid[i * _gridSize.y + _gridSize.y-1].soundPressure;
	}
	for(int j=0;j<_gridSize.y;j++){
		_murX1[0 + j * 4] = _grid[0 * _gridSize.y + j].soundPressure;
		_murX1[1 + j * 4] = _grid[1 * _gridSize.y + j].soundPressure;
		_murX1[2 + j * 4] = _grid[(_gridSize.x-2) * _gridSize.y + j].soundPressure;
		_murX1[3 + j * 4] = _grid[(_gridSize.x-1) * _gridSize.y + j].soundPressure;
	}
}

void AcousticFDTD::draw()
{
	updateV();
	updateP();

	mur2nd();

	/* Initial Waveform from a Point Source (1 pulse of sinusoidal wave with Hann window) */
	if( _nPoint < (1.0/_freq)/_dt ){
		_sigPoint = (1.0-cos((2.0*M_PI*_freq*_nPoint*_dt)))/2.0 * sin((2.0*M_PI*_freq*_nPoint*_dt));
		_grid[_pointSource.x * _gridSize.y + _pointSource.y].soundPressure += _sigPoint;
	}

	//setGridColors
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

	++_nPoint;

	if(_nPoint > _randomPointSourceInterval)
	{
		randomPointSource();
	}
}

void AcousticFDTD::randomPointSource()
{
	_pointSource.x = rand() % _gridSize.x;
	_pointSource.y = rand() % _gridSize.y;
	_nPoint = 0;
}
