#pragma once

#include "Tomo_types.h"
#include "STomoPixel.h"

using namespace Tomo;

template <class T> TOMO_LONGINT  _TomoNV_Function_Call(TOMO_FLOAT32* _info, TOMO_FLOAT32* _YPR, TOMO_FLOAT32* chull_vtx, TOMO_FLOAT32* _vtx, TOMO_FLOAT32* _tri_nrm, TOMO_LONGINT* _tri);
template <typename T> void thread_func(T* _pNV, int thread_id, TOMO_FLOAT32* _YPR, int ypr_id);//멀티스레드 적용을 위한 함수

TOMO_LONGINT _find1stOptimal(TOMO_LONGINT _nData, TOMO_FLOAT32* _pData);

#ifdef __cplusplus
extern "C"
{
#endif

DLLEXPORT TOMO_INT16*       getpData2i(Tomo::enumPixelType iPixel);
DLLEXPORT TOMO_LONGINT  getnData2i(Tomo::enumPixelType iPixel);
DLLEXPORT TOMO_FLOAT32* getMss(void);
DLLEXPORT TOMO_FLOAT32* getMo(void);
DLLEXPORT TOMO_FLOAT32* getVtc(void);//for p-orbital. python version.
DLLEXPORT TOMO_FLOAT32* getVolMassInfo(void);

DLLEXPORT TOMO_INT16* pxlsToDat2i(TPVector& pxls, TOMO_LONGINT& n_pxl);

DLLEXPORT void  OnDestroy(void);

//python interface functions
DLLEXPORT TOMO_LONGINT  TomoNV_TMPxl(TOMO_FLOAT32* _info, TOMO_FLOAT32* _YPR, TOMO_FLOAT32* chull_vtx, TOMO_FLOAT32* _vtx, TOMO_FLOAT32* _tri_nrm, TOMO_LONGINT* _tri);

#ifdef __cplusplus
}
#endif
