#pragma once
#include "STomoPixel.h"
#include <functional> //std::function

using namespace Tomo;
class DLLEXPORT STPSlot
{
public:
  STPSlot();
  STPSlot(const STPSlot& Source);
    void	operator=(const STPSlot& Source);
    void	_Copy(const STPSlot& Source);
  ~STPSlot();

  void	Reset(void);
  void	Init(void);

  TPVector  pxls;//sums of the following *_pxls.
    TPVector al_pxls, be_pxls, TC_pxls, NVB_pxls, NVA_pxls, Vo_pxls, Vss_pxls, SS_pxls, SSB_pxls, SSA_pxls;//internal variables for calculation.

  static const int niData = 2;
  static const int siData = sizeof(TOMO_INT16) * niData;
  union
  { struct  {  TOMO_INT16 X, Y;};//slot 2D coordinates.
                TOMO_INT16 iData[niData];
  };

  STomoVolMassInfo  vm_info;//local volume/mass values
  std::function<TOMO_FLOAT32(TOMO_FLOAT32)> wFunc;//obsolete..

  void  Pairing(TOMO_FLOAT32 theta_c, bool _bUseExplicitSS = false);//chapter 3.3 & 3.4
    TPVector  createTCPxls(void);
    void  createAlBePxls(void);
    TPVector createVoPxls(const TPVector& alpha, const TPVector& beta);
    void  splitAlBe(const TPVector& pxls, TPVector& al_pxls, TPVector& be_pxls);
  void  Calculate(void);

  //for rendering
  TPVector  slotToPxls(enumPixelType _type);
  TPVector  GetSSPxls(bool _bUseExplicitSS = false); 
};

  
DLLEXPORT typedef std::vector<STPSlot>		      TPSlotVector;
DLLEXPORT typedef TPSlotVector::iterator				TPSlotIterator;
DLLEXPORT typedef TPSlotVector::const_iterator	TPSlotConstIterator;

