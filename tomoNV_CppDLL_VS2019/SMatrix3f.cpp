#include "pch.h"
#include "SMatrix3f.h"
#include <cstring> //memcpy()
#include <math.h>

using namespace Tomo;

namespace Tomo
{

	SMatrix3f::SMatrix3f(_FLOAT yaw, _FLOAT pitch, _FLOAT roll)
	{
		YPR(yaw, pitch, roll);
	}

	SMatrix3f::~SMatrix3f()
	{
		Reset();
	}

	SMatrix3f::SMatrix3f(const SMatrix3f& Source)
	{
		Init();
		_Copy(Source);
	}

	void	SMatrix3f::operator=(const SMatrix3f& Source)
	{
		Reset();
		_Copy(Source);
	}

	void	SMatrix3f::_Copy(const SMatrix3f& Source)
	{
		memcpy(dData, Source.dData, sdData);
	}

	void	SMatrix3f::Reset(void)
	{
		Init();
	}

	void	SMatrix3f::Init(void)
	{
		memset(dData, 0x00, sdData);
	}

	//-------------------------------------------------
	SMatrix3f& SMatrix3f::Identity(void)
	{
		Reset();
		Data[0][0] = 1.;
		Data[1][1] = 1.;
		Data[2][2] = 1.;

		return *this;
	}

	SMatrix3f& SMatrix3f::OuterProduct( _FLOAT *a, _FLOAT *b)
	{
		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				Data[i][j] = a[i] * b[j];
			}
		}

		return *this;
	}

	SMatrix3f	SMatrix3f::operator+(const SMatrix3f&	a) const
	{
		SMatrix3f temp;

		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				temp.Data[i][j] = Data[i][j] + a.Data[i][j];
			}
		}

		return temp;
	}

	SMatrix3f	SMatrix3f::operator-(const SMatrix3f&	a) const
	{
		SMatrix3f temp;

		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				temp.Data[i][j] = Data[i][j] - a.Data[i][j];
			}
		}

		return temp;
	}

	SMatrix3f& SMatrix3f::operator+=(const SMatrix3f& a)
	{
		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				Data[i][j] += Data[i][j];
			}
		}

		return *this;
	}

	SMatrix3f& SMatrix3f::operator-=(const SMatrix3f& a)
	{
		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				Data[i][j] -= a.Data[i][j];
			}
		}

		return *this;
	}

	SMatrix3f SMatrix3f::operator*(_FLOAT b)
	{
		SMatrix3f temp;

		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				temp.Data[i][j] = Data[i][j] * b;
			}
		}

		return temp;
	}

	SMatrix3f operator*(_FLOAT b,SMatrix3f const & m)
	{
		SMatrix3f temp;

		for(int i = 0 ; i < 3 ; i++)
		{
			for( int j = 0 ; j < 3 ; j++)
			{
				temp.Data[i][j] = m.Data[i][j] * b;
			}
		}

		return temp;
	}

	void  SMatrix3f::Dot(_FLOAT* a, _FLOAT* result)
	{
		_FLOAT temp[3];

		temp[0]  = Data[0][0] * a[0];
		temp[0] += Data[0][1] * a[1];
		temp[0] += Data[0][2] * a[2];

		temp[1]  = Data[1][0] * a[0];
		temp[1] += Data[1][1] * a[1];
		temp[1] += Data[1][2] * a[2];

		temp[2]  = Data[2][0] * a[0];
		temp[2] += Data[2][1] * a[1];
		temp[2] += Data[2][2] * a[2];

		result[0] = temp[0];
		result[1] = temp[1];
		result[2] = temp[2];
	}


	void		SMatrix3f::T(SMatrix3f& a)
	{
		int i,j;
		for( i = 0 ; i < 3 ; i++)
		{
			for( j = 0 ; j < 3 ; j++)
			{
				Data[i][j] = a.Data[j][i];
			}
		}
	}

	void		SMatrix3f::AXPY( _FLOAT alpha, const SMatrix3f& a, _FLOAT a_factor)
	{
		int i,j;
		for( i = 0 ; i < 3 ; i++)
		{
			for( j = 0 ; j < 3 ; j++)
			{
				Data[i][j] = Data[i][j] * alpha + a.Data[i][j] * a_factor;
			}
		}
	}

	
	void		SMatrix3f::AXPY( _FLOAT a, float *x, _FLOAT b)
	{//for UpperDiagonalJacobian
		int i,j;
		for( i = 0 ; i < 3 ; i++)
		{
			for( j = 0 ; j < 3 ; j++)
			{
				Data[i][j] += a * x[i*3+j] + b;
			}
		}

	}

	void	  SMatrix3f::AddProduct( _FLOAT alpha, _FLOAT beta, const SMatrix3f& a, const SMatrix3f& b)
	{
		int i,j,k;
		_FLOAT sum;
		for( i = 0 ; i < 3 ; i++)
		{
			sum = 0;
			for( j = 0 ; j < 3 ; j++)
			{
				for( k = 0 ; k< 3 ; k++)
				{
					sum += a.Data[i][k] * b.Data[k][j];
				}
				Data[i][j] = alpha * Data[i][j] + beta * sum;
			}
		}
	}

	void    SMatrix3f::YPR(_FLOAT ga, _FLOAT be, _FLOAT al)//switch al, ga to match with python version. 
	{//http://planning.cs.uiuc.edu/node102.html
	 
		_FLOAT cos_al = ::cos(al);
		_FLOAT cos_be = ::cos(be);
		_FLOAT cos_ga = ::cos(ga);

		_FLOAT sin_al = ::sin(al);
		_FLOAT sin_be = ::sin(be);
		_FLOAT sin_ga = ::sin(ga);

		Data[0][0] = cos_al * cos_be;
		Data[0][1] = cos_al * sin_be * sin_ga - sin_al*  cos_ga;
		Data[0][2] = cos_al * sin_be * cos_ga + sin_al * sin_ga;

		Data[1][0] = sin_al * cos_be;
		Data[1][1] = sin_al * sin_be * sin_ga + cos_al * cos_ga;
		Data[1][2] = sin_al * sin_be * cos_ga - cos_al * sin_ga;

		Data[2][0] = - sin_be;
		Data[2][1] = cos_be * sin_ga;
		Data[2][2] = cos_be * cos_ga;
	}

}