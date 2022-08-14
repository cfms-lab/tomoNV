#pragma once
#include "Tomo_types.h"

class DLLEXPORT SClass
{
public:
  SClass();
  SClass(const SClass& Source);
  void	operator=(const SClass& Source);
  void	_Copy(const SClass& Source);
  ~SClass();

  void	Reset(void);
  void	Init(void);

  static const int niData = 6;
  static const int siData = sizeof(int) * niData;
  union
  {
    struct { int	x, y, z, nx, ny, nz; };
    int iData[niData];
  };

};

