// TomoNVC_Win32.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "TomoNV_Win64.h"
#include "STomoPixel.h"
#include "STPSlot.h"
#include "STomoNV_TMPxl.h"

#include <algorithm>  // std::find_if
#include <thread>
#include <iostream>//std::cout

#include "SMatrix3f.h"

//global variables
TOMO_FLOAT32** vtx = nullptr;
TOMO_FLOAT32** nrm = nullptr;
TOMO_INT16  ** tri = nullptr;
TOMO_INT16 nV = 0;
TOMO_INT16 nT = 0;
TOMO_INT16 n_pxls = 0;
TOMO_LONGINT*  nData2i = nullptr;
TOMO_INT16** pData2i = nullptr;
TOMO_FLOAT32* Mss = nullptr;
TOMO_FLOAT32* Mo = nullptr;
TOMO_FLOAT32* Vtc = nullptr;//debug. for p-orbital test
TOMO_FLOAT32 optimal_YPR[3];
bool  g_bUseExplicitSS = false;
STomoVolMassInfo VolMassInfo;

TPVector al_pxls, be_pxls, TC_pxls, NVB_pxls, NVA_pxls, Vo_pxls, Vss_pxls, SS_pxls, SSB_pxls, SSA_pxls;

using namespace Tomo;

inline TOMO_FLOAT32 toRadian(TOMO_INT16 a) { return a * 3.141592 / 180.; }

void triPixel(
  TOMO_FLOAT32* v0, TOMO_FLOAT32* v1, TOMO_FLOAT32* v2,
  TOMO_FLOAT32* n0, TOMO_FLOAT32* n1, TOMO_FLOAT32* n2,
  TPVector& tri_pxls)
{
  using namespace Tomo;
 
  int x0 = int(Tomo::_min(Tomo::_min(v0[0], v1[0]), v2[0]));
  int y0 = int(Tomo::_min(Tomo::_min(v0[1], v1[1]), v2[1]));
  int x1 = int(Tomo::_max(Tomo::_max(v0[0], v1[0]), v2[0]));
  int y1 = int(Tomo::_max(_max(v0[1], v1[1]), v2[1]));

  TOMO_FLOAT32 v_center[3] = { 0.,0.,0. }, u, v, w;

  giveMargin(v0);
  giveMargin(v1);
  giveMargin(v2);

  for (int x = x0; x <= x1 + 1; x++)
  {
    v_center[0] = TOMO_FLOAT32(x + g_fMARGIN);
    for (int y = y0; y <= y1 + 1; y++)
    {
      v_center[1] = TOMO_FLOAT32(y + g_fMARGIN);
      if (_getBaryCoord(v_center, v0, v1, v2, u, v, w))
      {
        TOMO_FLOAT32 pxl[3], nrm[3];
        _bary_product(v0, v1, v2, u, v, w, pxl);
        _bary_product(n0, n1, n2, u, v, w, nrm);
        STomoPixel new_pxl(pxl, nrm);
        tri_pxls.push_back(new_pxl);
      }
    }
  }
}



TOMO_INT16* pxlsToDat2i(TPVector& pxls, TOMO_LONGINT& n_pxl)
{
  n_pxl = pxls.size();
  if(n_pxl<=0) return nullptr;

  TOMO_INT16* _Data2i = new TOMO_INT16[n_pxl * g_nPixelFormat];
  memset(_Data2i, 0x00, sizeof(TOMO_INT16) * n_pxl * g_nPixelFormat);
  TOMO_LONGINT p = 0;
  for (TPIterator pIt = pxls.begin(); pIt != pxls.end(); ++pIt, ++p)
  {
    pIt->DumpTo(_Data2i + p * g_nPixelFormat);
  }
  return _Data2i;
}

TOMO_INT16*       getpData2i(Tomo::enumPixelType iSubPixel) { return ::pData2i[static_cast<int>(iSubPixel)];}
TOMO_LONGINT  getnData2i(Tomo::enumPixelType iSubPixel) { return ::nData2i[static_cast<int>(iSubPixel)];}
TOMO_FLOAT32* getMss(void) { return ::Mss;}
TOMO_FLOAT32* getMo(void)  { return ::Mo; }
TOMO_FLOAT32* getVtc(void) { return ::Vtc; }
TOMO_FLOAT32* getVolMassInfo(void) 
{ 
  return ::VolMassInfo.dData; 
}

void  OnDestroy(void)
{
  if(pData2i!=nullptr)
  {
    for (int i = 0; i < static_cast<int>(enumPixelType::espNumberOfSubPixels); i++)    {      if(pData2i[i] != nullptr) delete[] pData2i[i];    }
    delete[] pData2i;    pData2i = nullptr;
  }
  vtx = nullptr;  nrm = nullptr;  tri = nullptr;
  if (Mss != nullptr) { delete[] Mss;  Mss = nullptr;  }
  if (Mo != nullptr)  { delete[] Mo;   Mo = nullptr; }
  if (Vtc != nullptr) { delete[] Vtc;  Vtc = nullptr; }
}

TOMO_LONGINT _find1stOptimal(TOMO_LONGINT _nData, TOMO_FLOAT32* _pData)
{
  TOMO_LONGINT min_index = 0;
  TOMO_FLOAT32 min_value = TOMO_FLOAT32( 1e5);
  for (TOMO_LONGINT i = 0; i < _nData; i++)  {  if (_pData[i] < min_value)  { min_value = _pData[i]; min_index = i;  }  }
  return min_index;
}

template <typename T> void thread_func(T* _pNV, int thread_id, TOMO_FLOAT32* _YPR, int ypr_id)
{
  T* nv = _pNV + thread_id;

  S3DPrinterInfo& P_info = nv->printer_info;
  P_info.yaw = _YPR[ypr_id * 3 + 0];//input YPR data here.
  P_info.pitch = _YPR[ypr_id * 3 + 1];
  P_info.roll = _YPR[ypr_id * 3 + 2];

  nv->Rotate();
  nv->Pixelize();
  nv->Pairing();
  nv->Calculate();

  Mo[ypr_id] = nv->vm_info.Mo;
  Mss[ypr_id] = nv->vm_info.Mss;
  Vtc[ypr_id] = nv->vm_info.Vtc;
}

template <class T>
TOMO_LONGINT  _TomoNV_Function_Call(TOMO_FLOAT32* _info, TOMO_FLOAT32* _YPR, TOMO_FLOAT32* _chull_vtx, 
TOMO_FLOAT32* _vtx, TOMO_FLOAT32* _tri_nrm, TOMO_LONGINT* _tri)
{
  bool  bVerbose = (_info[0] > 0.);
  bool  bUseExplicitSS = (_info[1] > 0.);
  TOMO_FLOAT32 dVoxel = _info[2];//voxel size.1 mm
  int   nVoxel = int(_info[3]);//numver of voxels. 256^3
  TOMO_LONGINT nYPR = TOMO_LONGINT(_info[4]);//number of orientations to search. (y,p,r) / 3

  Mss = new TOMO_FLOAT32[nYPR + 2];
  Mo = new TOMO_FLOAT32[nYPR + 2];
  Vtc = new TOMO_FLOAT32[nYPR + 2];

  using namespace Tomo;

  //mult-thread info.
  const auto processor_count = std::thread::hardware_concurrency();
  int nThread = min(processor_count - 1, nYPR);
  int nBlock = nYPR / nThread;
  int nBlRest = nYPR % nThread;

  S3DPrinterInfo info;
  info.Set(_info, _vtx, _tri_nrm, _tri, 0, 0, 0);//takes some memory.

  T *pNV = new T[nThread +2];

  for (int thread_id = 0; thread_id < nThread; thread_id++)
  {
    pNV[thread_id].printer_info  = info;
  }

#ifdef _DEBUG
  //non-threaded version
  int ypr_id = 0;
  for (int b = 0; b < nBlock; b++)
  {
    for (int thread_id = 0; thread_id < nThread; thread_id++) { thread_func<T>(pNV, thread_id, _YPR, ypr_id++); }
    if (bVerbose && b % nBlock == 0) { std::cout << "Step" << ypr_id + 1 << "/" << nYPR << std::endl; }
  }

  { for (int thread_id = 0; thread_id < nBlRest; thread_id++) { thread_func<T>(pNV, thread_id, _YPR, ypr_id++); }  }
#else
  //std::thread version
  int ypr_id = 0;
  for (int b = 0; b < nBlock; b++) //Repeat # of CPU core * integer(nBlock) jobs
  {
    std::vector< std::thread> thVec;
    for (int thread_id = 0; thread_id < nThread; thread_id++) { thVec.push_back(std::thread(thread_func<T>, pNV, thread_id, _YPR, ypr_id++)); }
    for (auto& th : thVec) { th.join(); }

    if (bVerbose && b % 100 == 0) { std::cout << "Step" << ypr_id + 1 << "/" << nYPR << std::endl; }
}

  { std::vector< std::thread> thVec; //Do the rest jobs
  for (int thread_id = 0; thread_id < nBlRest; thread_id++) { thVec.push_back(std::thread(thread_func<T>, pNV, thread_id, _YPR, ypr_id++)); }
  for (auto& th : thVec) { th.join(); }  }

#endif

  TOMO_LONGINT  optID = (nYPR > 10) ? _find1stOptimal(nYPR, Mss) : 0;
  if(nYPR>1)  thread_func<T>(pNV, 0, _YPR, optID);//Find information of optID again.

  //prepare rendering data for python
  int nPixelType = static_cast<int>(enumPixelType::espNumberOfSubPixels);
  ::nData2i = new TOMO_LONGINT[nPixelType];
  ::pData2i = new TOMO_INT16*     [nPixelType];

  for (int i = 0; i < nPixelType; i++)
  {
    TPVector tmp_pxls = pNV[0].slotsToPxls(enumPixelType(i));
    pData2i[i] = pxlsToDat2i( tmp_pxls, nData2i[i]);
  }

  {//Find SS_pxls for rendering.
    TPVector tmp_pxls = pNV[0].GetSSPixels(bUseExplicitSS);
    int _SS = static_cast<int>(enumPixelType::eptSS);
    pData2i[_SS] = pxlsToDat2i(tmp_pxls, nData2i[_SS]);
  }

  VolMassInfo = pNV[0].vm_info; //final result

  delete[] pNV;
  return optID;
}

//python interface functions
TOMO_LONGINT  TomoNV_TMPxl(TOMO_FLOAT32* _info, TOMO_FLOAT32* _YPR, TOMO_FLOAT32* _chull_vtx, TOMO_FLOAT32* _vtx, TOMO_FLOAT32* _tri_nrm, TOMO_LONGINT* _tri)
{
  return  _TomoNV_Function_Call<STomoNV_TMPxl>(_info, _YPR, _chull_vtx, _vtx, _tri_nrm, _tri);
}
