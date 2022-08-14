#include "pch.h"
#include "Tomo_types.h"
#include "STomoPixel.h"
#include <vector>
#include <algorithm>  // std::find_if

using namespace Tomo;

namespace Tomo
{
  S3DPrinterInfo::S3DPrinterInfo()
  {
    Init();
  }

  S3DPrinterInfo::~S3DPrinterInfo()
  {
    Reset();
  }

  S3DPrinterInfo::S3DPrinterInfo(const S3DPrinterInfo& Source)
  {
    Init();
    _Copy(Source);
  }

  void	S3DPrinterInfo::operator=(const S3DPrinterInfo& Source)
  {
    Reset();
    _Copy(Source);
  }

  void	S3DPrinterInfo::_Copy(const S3DPrinterInfo& Source)
  {
    memcpy( dData, Source.dData, sdData);

    bVerbose = Source.bVerbose;
    bUseExplicitSS = Source.bUseExplicitSS;

    rpVtx0 = Source.rpVtx0;
    rpTriNrm0 = Source.rpTriNrm0;
    rpTri0 = Source.rpTri0;
    nVtx  = Source.nVtx;
    nTri  = Source.nTri;
    nVoxel = Source.nVoxel;
    nYPR  = Source.nYPR;

    //convexhull test
    nCHVertices = Source.nCHVertices;
    if(nCHVertices > 0 && Source.pCHVertices != nullptr)
    {
      pCHVertices = new TOMO_FLOAT32[nCHVertices*3+2];
      memcpy(pCHVertices, Source.pCHVertices, sizeof(TOMO_FLOAT32) * nCHVertices * 3);
    }

  }

  void	S3DPrinterInfo::Reset(void)
  {
    if (pVtx1  != nullptr) { delete[] pVtx1;   pVtx1 = nullptr; }
    if (pTriNrm1  != nullptr) { delete[] pTriNrm1;   pTriNrm1 = nullptr; }

    Init();
  }

  void	S3DPrinterInfo::Init(void)
  {
    memset(dData, 0x00, sdData);

    bVerbose = bUseExplicitSS = false;

    rpVtx0 = rpTriNrm0 = pVtx1 = pTriNrm1 = nullptr;
    rpTri0 = nullptr;
    nVtx = nTri = 0;
    nVoxel = 256;
    nYPR = 0;

    //convexhull test
    nCHVertices = 0;
    pCHVertices = nullptr;
  }

  void	S3DPrinterInfo::Set(TOMO_FLOAT32 *_info, 
    TOMO_FLOAT32* _vtx, TOMO_FLOAT32* _tri_nrm, TOMO_LONGINT* _tri,
    TOMO_FLOAT32 _yaw, TOMO_FLOAT32 _pitch, TOMO_FLOAT32 _roll )
  {
    bVerbose        = (_info[0] > g_fMARGIN);
    bUseExplicitSS  = (_info[1] > g_fMARGIN);
    dVoxel  = _info[2];
    nVoxel  = int(_info[3]);
    nYPR  = TOMO_INT16(_info[4]);

    theta_c = _info[5];

    nVtx    = TOMO_LONGINT(_info[6]);
    nTri    = TOMO_LONGINT(_info[7]);

    surface_area    = _info[8];
    wall_thickness  = _info[9];
    PLA_density     = _info[10];

    Fcore = _info[11];
    Fclad = _info[12];
    Fss   = _info[13];
    Css   = _info[14];

    yaw   = _yaw;
    pitch = _pitch;
    roll  = _roll;

    rpVtx0  = _vtx;
    rpTriNrm0  = _tri_nrm;
    pVtx1   = nullptr;
    rpTri0  = _tri;

    nCHVertices = int(_info[15]);
  }



  void  STomoAABB2D::GetAABB2D(TOMO_LONGINT _nVtx, TOMO_FLOAT32* _pVtx)
  {
    x0 = TOMO_INT16(1e4);  x1 = TOMO_INT16(-1e4);
    y0 = TOMO_INT16(1e4);  y1 = TOMO_INT16(-1e4);
    for (TOMO_LONGINT v = 0; v < _nVtx; v++)
    {
      TOMO_INT16 x = TOMO_INT16( _pVtx[v * 3 + 0]); //FLOAT64 -> INT64
      TOMO_INT16 y = TOMO_INT16( _pVtx[v * 3 + 1]); //FLOAT64 -> INT64
      x0 = _min(x, x0);
      y0 = _min(y, y0);
      x1 = _max(x, x1);
      y1 = _max(y, y1);
    }
  }

  //------------------------------------------------------------------

  TOMO_FLOAT32 inline _dot(TOMO_FLOAT32* a, TOMO_FLOAT32* b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  TOMO_FLOAT32 inline _abs(TOMO_FLOAT32 value)
  {
    return (value >= 0.) ? value : -value;
  }

  TOMO_INT16 inline _abs(TOMO_INT16 value)
  {
    return (value >= 0) ? value : -value;
  }

  TOMO_FLOAT32 inline _min(TOMO_FLOAT32 a, TOMO_FLOAT32 b)
  {
    return (a <= b) ? a : b;
  }

  TOMO_INT16 inline _min(TOMO_INT16 a, TOMO_INT16 b)
  {
    return (a <= b) ? a : b;
  }

  TOMO_FLOAT32 inline _max(TOMO_FLOAT32 a, TOMO_FLOAT32 b)
  {
    return (a >= b) ? a : b;
  }

  TOMO_INT16 inline _max(TOMO_INT16 a, TOMO_INT16 b)
  {
    return (a >= b) ? a : b;
  }

  SLOT_BUFFER_TYPE   inline _round(TOMO_FLOAT32 a)
  {
    int a10 = int(a * 10);
    int digit1 = a10 - int(a)*10;
    if(digit1>=5) return int(a+1);
    else          return int(a);
  }

  DLLEXPORT TOMO_FLOAT32 inline _toDegree(TOMO_FLOAT32 _radian)
  {
    return TOMO_FLOAT32(_radian * 180. / 3.141592);
  }

  DLLEXPORT TOMO_FLOAT32 inline _toRadian(TOMO_FLOAT32 _degree)
  {
    return TOMO_FLOAT32(_degree / 180. * 3.141592);
  }


  void inline giveMargin(TOMO_FLOAT32* v)
  {
    v[0] += g_fMARGIN;
    v[1] += g_fMARGIN;
    v[2] += g_fMARGIN;
  }

  void inline _bary_product(
    /*inputs*/ TOMO_FLOAT32* p0, TOMO_FLOAT32* p1, TOMO_FLOAT32* p2,
    TOMO_FLOAT32  u, TOMO_FLOAT32   v, TOMO_FLOAT32   w,
    /*output*/ TOMO_FLOAT32* pxl)
  {
    for (int i = 0; i < 3; i++)
    {
      pxl[i] = p0[i] * u + p1[i] * v + p2[i] * w;
    }
  }

  bool  _getBaryCoord(
    /*inputs*/ TOMO_FLOAT32* p, TOMO_FLOAT32* a, TOMO_FLOAT32* b, TOMO_FLOAT32* c,
    /*output*/ TOMO_FLOAT32& u, TOMO_FLOAT32& v, TOMO_FLOAT32& w)
  {
    TOMO_FLOAT32 v0[3] = { b[0] - a[0], b[1] - a[1], 0. };
    TOMO_FLOAT32 v1[3] = { c[0] - a[0], c[1] - a[1], 0. };
    TOMO_FLOAT32 v2[3] = { p[0] - a[0], p[1] - a[1], 0. };

    TOMO_FLOAT32 d00 = _dot(v0, v0);
    TOMO_FLOAT32 d01 = _dot(v0, v1);
    TOMO_FLOAT32 d11 = _dot(v1, v1);
    TOMO_FLOAT32 d20 = _dot(v2, v0);
    TOMO_FLOAT32 d21 = _dot(v2, v1);

    TOMO_FLOAT32 denom = d00 * d11 - d01 * d01;

    if (_abs(denom) > g_fMARGIN)
    {
      v = (d11 * d20 - d01 * d21) / denom;
      w = (d00 * d21 - d01 * d20) / denom;
      u = 1.0 - v - w;
      if (u >= -g_fMARGIN && v >= -g_fMARGIN && v <= 1. + g_fMARGIN && u + v <= 1. + g_fMARGIN)
      {
        return true;
      }
    }
    //else
    u = -1.; v = -1.; w = -1.;
    return false;
  }


  void  inline QuickSort(SLOT_BUFFER_TYPE A[], size_t I[], size_t lo, size_t hi)
  {//thanks to https://stackoverflow.com/questions/55976487/get-the-sorted-indices-of-an-array-using-quicksort
    while (lo < hi)
    {
      SLOT_BUFFER_TYPE pivot = A[I[lo + (hi - lo) / 2]];
      size_t t;
      size_t i = lo - 1;
      size_t j = hi + 1;
      while (1)
      {
        while (A[I[++i]] < pivot);
        while (A[I[--j]] > pivot);
        if (i >= j)
          break;
        t = I[i];
        I[i] = I[j];
        I[j] = t;
      }
      /* avoid stack overflow */
      if ((j - lo) < (hi - j)) {
        QuickSort(A, I, lo, j);
        lo = j + 1;
      }
      else {
        QuickSort(A, I, j + 1, hi);
        hi = j;
      }
    }
  }

  SLOT_BUFFER_TYPE  inline toTypeByte(enumPixelType type)
  {
    return SLOT_BUFFER_TYPE(1 << static_cast<SLOT_BUFFER_TYPE>(type));
  }

  //------------------------------------------------------------------

  STomoAABB3Df::STomoAABB3Df()
  {
    Init();
  }

  STomoAABB3Df::~STomoAABB3Df()
  {
    Reset();
  }

  STomoAABB3Df::STomoAABB3Df(const STomoAABB3Df& Source)
  {
    Init();
    _Copy(Source);
  }

  void	STomoAABB3Df::operator=(const STomoAABB3Df& Source)
  {
    Reset();
    _Copy(Source);
  }

  void	STomoAABB3Df::_Copy(const STomoAABB3Df& Source)
  {
    memcpy(dData, Source.dData, sdData);
  }

  void	STomoAABB3Df::Reset(void)
  {
    Init();
  }

  void	STomoAABB3Df::Init(void)
  {
    //memset(dData, 0x00, sdData);
    x0 = y0 = z0 = TOMO_FLOAT32(1e5);
    x1 = y1 = z1 = TOMO_FLOAT32(-1e5);
  }

  void  STomoAABB3Df::Set(TOMO_LONGINT nV, TOMO_FLOAT32* vtx)
  {
    Init();
    for (TOMO_LONGINT v = 0; v < nV; v++)
    {
      TOMO_FLOAT32 x = *(vtx + v * 3 + 0);
      TOMO_FLOAT32 y = *(vtx + v * 3 + 1);
      TOMO_FLOAT32 z = *(vtx + v * 3 + 2);
      x0 = _min(x, x0);
      y0 = _min(y, y0);
      z0 = _min(z, z0);

      x1 = _max(x, x1);
      y1 = _max(y, y1);
      z1 = _max(z, z1);
    }
  }

  void  STomoAABB3Df::GetCenter(TOMO_FLOAT32* _center)
  {
    _center[0] = TOMO_FLOAT32((x0 + x1) * 0.5);
    _center[1] = TOMO_FLOAT32((y0 + z1) * 0.5);
    _center[2] = TOMO_FLOAT32((y0 + z1) * 0.5);
  }
  //------------------------------------------------------------------
  STomoVolMassInfo::STomoVolMassInfo()
  {
    Init();
  }

  STomoVolMassInfo::~STomoVolMassInfo()
  {
    Reset();
  }

  STomoVolMassInfo::STomoVolMassInfo(const STomoVolMassInfo& Source)
  {
    Init();
    _Copy(Source);
  }

  void	STomoVolMassInfo::operator=(const STomoVolMassInfo& Source)
  {
    Reset();
    _Copy(Source);
  }

  void	STomoVolMassInfo::_Copy(const STomoVolMassInfo& Source)
  {
    memcpy(dData, Source.dData, sdData);
  }

  void	STomoVolMassInfo::Reset(void)
  {
    Init();
  }

  void	STomoVolMassInfo::Init(void)
  {
    memset(dData, 0x00, sdData);
  }

}

void  STomoVolMassInfo::VolToMass(const S3DPrinterInfo&printer_info) //Eq. (10) ~ (16)
{
  Vclad = printer_info.surface_area * printer_info.wall_thickness;
  Vcore = Vo - Vclad;
  Mclad = Vclad * printer_info.Fclad * printer_info.PLA_density;
  Mcore = Vcore * printer_info.Fcore * printer_info.PLA_density;

#ifdef _EZAIR_THETAC_ZERO
  Vss = Vtc - Vo;
#endif

  Mo = Mcore + Mclad;
  Mss = Vss * printer_info.Fss * printer_info.Css * printer_info.PLA_density;
  Mtotal = Mo + Mss;
}
