#pragma once
#pragma warning ( disable:4251)
#pragma warning ( disable:4819)

#ifdef _CREATING_DLL_
#define DLLEXPORT __declspec(dllexport)
#define DLLEXTERN
#else
#define TomoNVCDLL_EXPORT __declspec(dllimport)
#define DLLEXTERN extern
#endif

//#define _EZAIR_THETAC_ZERO //Vss = Vtc - Vo. Use zero critical angle for test folliwing Ezair et al's paper.
#define _USE_VTX_NRM_FOR_PIXEL //user per-pixel Gouraud normal vector instead of flat normal. Decreases Vss value error slightly.

namespace Tomo 
{

  DLLEXPORT typedef float       TOMO_FLOAT32;//for mesh coordinates
  DLLEXPORT typedef short int   TOMO_INT16;
  DLLEXPORT typedef long int    TOMO_LONGINT;//for mesh element index

  DLLEXPORT typedef TOMO_INT16       SLOT_BUFFER_TYPE;//for slot internal data

  DLLEXPORT const int g_nPixelFormat = 6;// size of pixel data:  (x,y,z, nx, ny, nz)
  DLLEXPORT const TOMO_FLOAT32 g_fNORMALFACTOR = TOMO_FLOAT32(1000.); //multiply normal-z value to reduce floating-point-to-integer conversion error
  DLLEXPORT const TOMO_FLOAT32 g_fMARGIN = TOMO_FLOAT32(0.001);//Machine epsilon for floating point comparison 

  //DLLEXPORT const int iVOXELFACTOR = 1;//1==256 voxel. 2== 512 voxel... //Table 4. Testing the effect of voxel size. 

  DLLEXPORT enum class enumPixelType { //type of each pixels
    eptAl   = 0, 
    eptBe   = 1, 
    eptTC   = 2, 
    eptNVB  = 3, 
    eptNVA  = 4, 
    eptVo   = 5, 
    eptVss  = 6,
    eptSS   = 7,
    eptSSB  = 8,
    eptSSA  = 9,
    espNumberOfSubPixels};

const SLOT_BUFFER_TYPE typeAl = 1 << (int)enumPixelType::eptAl;//1
const SLOT_BUFFER_TYPE typeBe = 1 << (int)enumPixelType::eptBe;//2
const SLOT_BUFFER_TYPE typeTC = 1 << (int)enumPixelType::eptTC;//4
const SLOT_BUFFER_TYPE typeNVB = 1 << (int)enumPixelType::eptNVB;//8
const SLOT_BUFFER_TYPE typeNVA = 1 << (int)enumPixelType::eptNVA;//16
const SLOT_BUFFER_TYPE typeVo = 1 << (int)enumPixelType::eptVo;//32
const SLOT_BUFFER_TYPE typeVss = 1 << (int)enumPixelType::eptVss;//64
const SLOT_BUFFER_TYPE typeSS = 1 << (int)enumPixelType::eptSS;//128
const SLOT_BUFFER_TYPE typeSSB = 1 << (int)enumPixelType::eptSSB;//256
const SLOT_BUFFER_TYPE typeSSA = 1 << (int)enumPixelType::eptSSA;//512


  class DLLEXPORT STomoAABB2D //find the 2D boundary of pixels
  {
  public:
    STomoAABB2D(TOMO_INT16 _x0=0, TOMO_INT16 _y0=0, TOMO_INT16 _x1 = 0, TOMO_INT16 _y1 = 0)
          : x0(_x0), y0(_y0), x1(_x1), y1(_y1) {}
      STomoAABB2D(TOMO_INT16* aabb2d) :
      x_min(aabb2d[0]),
      y_min(aabb2d[1]),
      x_max(aabb2d[2]),
      y_max(aabb2d[3]) {}
    ~STomoAABB2D() {}

    static const int niData = 4;
    static const int siData = sizeof(TOMO_INT16) * niData;
    union
    {
      struct { TOMO_INT16	x_min, y_min, x_max, y_max; };
      struct { TOMO_INT16	x0, y0, x1, y1; };
      TOMO_INT16 iData[niData];
    };

    TOMO_INT16 nRow(void) const { return x1 - x0 + 1; }
    TOMO_INT16 nCol(void) const { return y1 - y0 + 1; }

    void  GetAABB2D(TOMO_LONGINT _nVtx, TOMO_FLOAT32* _pVtx);
  };

  class DLLEXPORT STomoAABB3Df //needed for mesh rotation operation.
  {
  public:
    STomoAABB3Df();
    STomoAABB3Df(const STomoAABB3Df& Source);
    void	operator=(const STomoAABB3Df& Source);
    void	_Copy(const STomoAABB3Df& Source);
    ~STomoAABB3Df();

    void	Reset(void);
    void	Init(void);

    static const int ndData = 6;
    static const int sdData = sizeof(TOMO_FLOAT32) * ndData;
    union
    {
      struct { TOMO_FLOAT32	x_min, y_min, x_max, y_max, z_min, z_max; };
      struct { TOMO_FLOAT32	x0, y0, x1, y1, z0, z1; };
      TOMO_FLOAT32 dData[ndData];
    };

    void  Set(TOMO_LONGINT nV, TOMO_FLOAT32* vtx);
    void  GetCenter(TOMO_FLOAT32 *center);
  };

  class DLLEXPORT S3DPrinterInfo //contains g-code parameters
  {
  public:
    S3DPrinterInfo();
    S3DPrinterInfo(const S3DPrinterInfo& Source);
      void	operator=(const S3DPrinterInfo& Source);
      void	_Copy(const S3DPrinterInfo& Source);
    ~S3DPrinterInfo();

    static const int ndData = 12;
    static const int sdData = sizeof(TOMO_FLOAT32) * ndData;
    union
    {
      struct { 
        TOMO_FLOAT32	dVoxel, theta_c, yaw, pitch, roll, \
          surface_area, wall_thickness, PLA_density, \
          Fcore, Fclad, Fss, Css; };
      TOMO_FLOAT32 dData[ndData];
    };

    bool  bVerbose;//for debug
    bool  bUseExplicitSS;

    //input mesh data
    TOMO_FLOAT32 *rpVtx0, *pVtx1;//[nVtx]
    TOMO_FLOAT32 *rpTriNrm0, * pTriNrm1;//[nTri]  triangle normal before/after rotation.
    TOMO_LONGINT *rpTri0;//[nTri]
    
    TOMO_LONGINT nVtx, nTri;//Cuastion: Do not treat this as int.
    int  nVoxel;//default 256
    int  nYPR;

    //convexhull test
    int nCHVertices;
    TOMO_FLOAT32* pCHVertices;

    void	Reset(void);
    void	Init(void);
    void  Set(TOMO_FLOAT32* _info, TOMO_FLOAT32* _vtx, TOMO_FLOAT32* _nrm, TOMO_LONGINT* _tri,
      TOMO_FLOAT32 _yaw, TOMO_FLOAT32 _pitch, TOMO_FLOAT32 _roll  );
      
  };

  class DLLEXPORT STomoVolMassInfo //contains the final volume & mass values
  {
  public:
    STomoVolMassInfo();
    STomoVolMassInfo(const STomoVolMassInfo& Source);
    void	operator=(const STomoVolMassInfo& Source);
    void	_Copy(const STomoVolMassInfo& Source);
    ~STomoVolMassInfo();

    void	Reset(void);
    void	Init(void);

    void  VolToMass(const S3DPrinterInfo&printer_info);  //Eq. (10) ~ (16)

    static const int ndData = 16;
    static const int sdData = sizeof(TOMO_FLOAT32) * ndData;
    union
    {
      struct {
        TOMO_FLOAT32	Va, Vb, Vtc, Vnv, \
          Vss, Vo, Vclad, Vcore,  \
          Mclad, Mcore, Mo, Mss, \
          Mtotal,/*for debug*/SS_vol;
      };
      TOMO_FLOAT32 dData[ndData];
    };
  };

  DLLEXPORT TOMO_FLOAT32 inline _dot(TOMO_FLOAT32* a, TOMO_FLOAT32* b);
  DLLEXPORT TOMO_FLOAT32 inline _abs(TOMO_FLOAT32 value);
  DLLEXPORT TOMO_INT16   inline _abs(TOMO_INT16 value);
  DLLEXPORT TOMO_FLOAT32 inline _min(TOMO_FLOAT32 a, TOMO_FLOAT32 b);
  DLLEXPORT TOMO_FLOAT32 inline _max(TOMO_FLOAT32 a, TOMO_FLOAT32 b);
  DLLEXPORT TOMO_INT16   inline _min(TOMO_INT16 a, TOMO_INT16 b);
  DLLEXPORT TOMO_INT16   inline _max(TOMO_INT16 a, TOMO_INT16 b);
  DLLEXPORT SLOT_BUFFER_TYPE   inline _round(TOMO_FLOAT32 a);
  DLLEXPORT TOMO_FLOAT32 inline _toDegree(TOMO_FLOAT32 _radian);
  DLLEXPORT TOMO_FLOAT32 inline _toRadian(TOMO_FLOAT32 _degree);
  DLLEXPORT void inline giveMargin(TOMO_FLOAT32* v);
  DLLEXPORT void inline _bary_product(
    /*inputs*/ TOMO_FLOAT32* p0, TOMO_FLOAT32* p1, TOMO_FLOAT32* p2,
    TOMO_FLOAT32  u, TOMO_FLOAT32   v, TOMO_FLOAT32   w,
    /*output*/ TOMO_FLOAT32* pxl);
  DLLEXPORT bool  _getBaryCoord(
    /*inputs*/ TOMO_FLOAT32* p, TOMO_FLOAT32* a, TOMO_FLOAT32* b, TOMO_FLOAT32* c,
    /*output*/ TOMO_FLOAT32& u, TOMO_FLOAT32& v, TOMO_FLOAT32& w);
  DLLEXPORT void  inline QuickSort(SLOT_BUFFER_TYPE A[], size_t I[], size_t lo, size_t hi);

  DLLEXPORT TOMO_FLOAT32 wFunc_3DWOX_Def(TOMO_FLOAT32 x);//obsolete..
  DLLEXPORT SLOT_BUFFER_TYPE  inline toTypeByte(enumPixelType type);

 
}