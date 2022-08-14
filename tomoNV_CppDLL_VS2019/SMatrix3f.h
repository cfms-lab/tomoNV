#pragma once
#include "Tomo_types.h"

namespace Tomo
{
  typedef TOMO_FLOAT32 _FLOAT;

  class	DLLEXPORT SMatrix3f {
  public:
	  SMatrix3f(_FLOAT yaw = 0., _FLOAT pitch = 0., _FLOAT roll = 0.);
    SMatrix3f(const SMatrix3f& Source);
    void	operator=(const SMatrix3f& Source);
    void	_Copy(const SMatrix3f& Source);
    ~SMatrix3f();

    void	Reset(void);
    void	Init(void);


	  static const int nRow = 3;
	  static const int nCol = 3;

    static const int ndData = nRow*nCol;
    static const int sdData = sizeof(_FLOAT) * ndData;
    union
    {
      _FLOAT Data[nRow][nCol];
      _FLOAT dData[nRow*nCol];
    };
    ;
    SMatrix3f&	Identity(void);
    SMatrix3f&	OuterProduct( _FLOAT *a, _FLOAT *b);//3d vector

	  SMatrix3f	  operator+(const SMatrix3f&	a) const;
    SMatrix3f	  operator-(const SMatrix3f&	a) const;
    SMatrix3f&  operator+=(const SMatrix3f&);
    SMatrix3f&  operator-=(const SMatrix3f&);

	  SMatrix3f operator*(_FLOAT);
	  void  Dot(_FLOAT*a, _FLOAT*result);

    void		T(SMatrix3f& a);//transpose

    void		AXPY( _FLOAT alpha, const SMatrix3f& a, _FLOAT a_factor =1.);
    void		AXPY( _FLOAT a, _FLOAT *x, _FLOAT b);//UpperDiagonalJacobian¿ë
	  void	  AddProduct( _FLOAT alpha, _FLOAT beta, const SMatrix3f& a, const SMatrix3f& b);

    void    YPR( _FLOAT yaw, _FLOAT pitch, _FLOAT roll);

  };

  SMatrix3f operator*(_FLOAT,SMatrix3f const&);

}