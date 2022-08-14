#pragma once
#include <vector>//stl
#include "Tomo_types.h"
#include <functional>

using namespace Tomo;

class DLLEXPORT STomoPixel //contains each pixel's coordinates + normal vectors 
{ 
public:
  STomoPixel(TOMO_INT16 _x=0, TOMO_INT16 _y = 0, TOMO_INT16 _z = 0, 
      TOMO_INT16 _nx = 0, TOMO_INT16 _ny = 0, TOMO_INT16 _nz = 0) : 
      x(_x), y(_y), z(_z), nx(_nx), ny(_ny), nz(_nz), iTypeByte(0) { }
  STomoPixel(TOMO_FLOAT32* _pxl3d, TOMO_FLOAT32 *_nrm3d);
  STomoPixel(TOMO_INT16* _data_6i);
  STomoPixel(const STomoPixel& Source);
    void	operator=(const STomoPixel& Source);
    void	_Copy(const STomoPixel& Source);
    void	DumpTo(TOMO_INT16* _data_2i);
  ~STomoPixel();

  void	Reset(void);
  void	Init(void);

  static const int niData = 6;
  static const int siData = sizeof(TOMO_INT16) * niData;
  union
  { struct  {   TOMO_INT16	x, y, z, nx, ny, nz;   };//Eq. (7)
                TOMO_INT16 iData[niData];
  };

  SLOT_BUFFER_TYPE    iTypeByte;//type of this pixel : alpha, beta, SSB, SSA, and so on.
};


DLLEXTERN template class DLLEXPORT std::vector<STomoPixel>;
DLLEXPORT typedef std::vector<STomoPixel>		 TPVector;
DLLEXPORT typedef TPVector::reverse_iterator				 TPReverseIterator;
DLLEXPORT typedef TPVector::iterator				 TPIterator;
DLLEXPORT typedef TPVector::const_iterator	 TPConstIterator;
DLLEXPORT typedef TPVector**  TPVectorpp;

DLLEXPORT TPVector& operator<<(TPVector& lhs, const TPVector& rhs);
DLLEXPORT bool operator==(const STomoPixel&lhs, const STomoPixel& rhs);
DLLEXPORT TPIterator _find(TPVector& pxls, const STomoPixel& pixel);

//functions
TPVector createVoPixels(const TPVector& alpha, const TPVector& beta, const STomoAABB2D& AABB2D);
TPVector createShadowCastorPixels(TOMO_FLOAT32 theta_c_in_Radian, TPVector& _be_pxls, bool _bExplicitSS=false);
TPVector createNVAPixels(const TPVector& nvB0, const TPVector& alpha0, const STomoAABB2D& AABB2D);
  STomoPixel  getNVA(const STomoPixel& nvb, const TPVector& al_slot);
TPVector createTCPixels(const TPVector& _al_pxls, const STomoAABB2D& AABB2D);
TPVector createSSPixels(
    const TPVector& alpha, const TPVector& beta,
    const TPVector& tc0,
    const TPVector& nvB, const TPVector& nvA,
    const Tomo::STomoAABB2D& AABB2D);
  void inline sortSlotByZ(TPVector& slot, bool _Higher = false);

TPVector    removeZNearPxls(TPVector& pxls0);
  TPVector    _matchAlternationWRT(const TPVector& slot, const TPVector& ref_slot);
  void        _matchPairNumber(TPVector& al_slot, TPVector& be_slot);//slot은 z값에 대해 오름차순임.
  bool       _hasPxlBetween(const TPVector& slot, TOMO_INT16 z_low, TOMO_INT16 z_high);
  void  _matchPairWRT(TPVector& pxls, const TPVector& ref_pxls);

  bool  _checkHigherPixel(const STomoPixel& a, const STomoPixel& b);
  bool  _checkLowerPixel(const STomoPixel& a, const STomoPixel& b);
  bool  _checkUnique(const STomoPixel& a, const STomoPixel& b);
  bool  _checkZNear(const STomoPixel& a, const STomoPixel& b);
  bool  _checkBottomPixel(const STomoPixel& a);

  //pixel operation
  void  copyPxls(const TPVector& src, TPVector& target);

  void  negativeZ(TPVector& pxls);
  void  zero_nZ(TPVector& pxls);
  //void  sortSlotByZ(TPVector& slot, bool _Higher = false);
  void  uniqueIntList(TPVector& pxls0);
  void  deleteBottomPxls(TPVector& pxls);
  void  deleteNoise(TPVector& pxls, TOMO_INT16 z_min, TOMO_INT16 z_max);
  TPVector getHighestPxls(TPVectorpp slots0, const STomoAABB2D& AABB2D);
  TPVector addZs(TPVectorpp slots0, const STomoAABB2D& AABB2D,
    std::function<TOMO_FLOAT32(TOMO_FLOAT32)> wFunc = nullptr);
  DLLEXPORT STomoPixel _addZ(const TPVector& pxls, 
    std::function<TOMO_FLOAT32(TOMO_FLOAT32)> wFunc = nullptr);
  //DLLEXPORT STomoPixel _weighted_addZ(const TPVector& pxls, std::function<TOMO_FLOAT32(TOMO_FLOAT32)> wFunc);

  //slot operation
  TPVectorpp  newSlots(const STomoAABB2D& AABB2D);
  void              deleteSlots(TPVectorpp slots, const STomoAABB2D& AABB2D);
  TPVectorpp  pxlsToSlots(const TPVector& pxls, const STomoAABB2D& AABB2D);
  void              slotsToPxls(TPVectorpp slots, const STomoAABB2D& AABB2D, TPVector& pxls);
  void  removeZNearPxls(TPVectorpp slots, const STomoAABB2D& AABB2D);

  void  matchAlBeSlots(TPVectorpp& al_slots, TPVectorpp& be_slots, const STomoAABB2D& AABB2D);
  //TPVector  _matchAlternationWRT(const TPVector& slot, const TPVector& ref_slot);
  //bool      _hasPxlBetween(const TPVector& slot, TOMO_INT16 z_low, TOMO_INT16 z_high);
  //void      _matchPairNumber(TPVector& al_slot, TPVector& be_slot);
  void  matchSlotsWRT(TPVectorpp& al_slots, TPVectorpp& be_slots, const STomoAABB2D& AABB2D);
  //void        _matchPairWRT(TPVector& pxls, const TPVector& ref_pxls);

  void  getAABB3D(TOMO_INT16 nV, TOMO_FLOAT32 *vtx, STomoAABB3Df& AABB3D);
