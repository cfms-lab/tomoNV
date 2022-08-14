#pragma once
#include "STPSlot.h"
#include "STomoNV_Base.h"

using namespace Tomo;

class DLLEXPORT STomoNV_TMPxl : public STomoNV_Base
{
public:
  STomoNV_TMPxl();
  STomoNV_TMPxl(const STomoNV_TMPxl& Source);
    void	operator=(const STomoNV_TMPxl& Source);
    void	_Copy(const STomoNV_TMPxl& Source);
  ~STomoNV_TMPxl();

  void	Reset(void);
  void	Init(void);
  
  void  Pixelize(void);
    void triPixel(
    TOMO_FLOAT32* v0, TOMO_FLOAT32* v1, TOMO_FLOAT32* v2,
    TOMO_FLOAT32* n0, TOMO_FLOAT32* n1, TOMO_FLOAT32* n2,
    TPVector& tri_pxls);
    TPVector  slotsToPxls(enumPixelType _type);//for rendering. time consuming. 
    void      pxlsToSlots(TPVector& tri_pxls);
  void  Pairing(void);//slot paring.
  void  Calculate(void);//get Vss value

  TPVector GetSSPixels(bool _bUseExplicitSS);

};
