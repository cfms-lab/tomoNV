#include "pch.h"
#include "STomoNV_Base.h"
#include "SMatrix3f.h"

using namespace Tomo;

STomoNV_Base::STomoNV_Base()
{
  Init();
}

STomoNV_Base::~STomoNV_Base()
{
  Reset();
}

STomoNV_Base::STomoNV_Base(const STomoNV_Base& Source)
{
  Init();
  _Copy(Source);
}

void	STomoNV_Base::operator=(const STomoNV_Base& Source)
{
  Reset();
  _Copy(Source);
}

void	STomoNV_Base::_Copy(const STomoNV_Base& Source)
{
  printer_info = Source.printer_info;
  slotVec = Source.slotVec;
}

void	STomoNV_Base::Reset(void)
{
  slotVec.clear();
}

void	STomoNV_Base::Init(void)
{
  slotVec.clear();
}

void  STomoNV_Base::Rotate(void)
{
  if (printer_info.pVtx1 == nullptr) printer_info.pVtx1 = new TOMO_FLOAT32[printer_info.nVtx * 3 + 2];
  if (printer_info.pTriNrm1 == nullptr) 
  {
#ifdef _USE_VTX_NRM_FOR_PIXEL
    printer_info.pTriNrm1 = new TOMO_FLOAT32[printer_info.nVtx * 3 + 2];
#else
    printer_info.pTriNrm1 = new TOMO_FLOAT32[printer_info.nTri * 3 + 2];
#endif
  }
    
  TOMO_FLOAT32* V0 = printer_info.rpVtx0;//raw data before rotation. do not change these.
  TOMO_FLOAT32* N0 = printer_info.rpTriNrm0;//raw data before rotation. do not change these.
  TOMO_FLOAT32* Vtx = printer_info.pVtx1;//after rotation.
  TOMO_FLOAT32* Nrm = printer_info.pTriNrm1;//after rotation.

  STomoAABB3Df AABB3Df;
  TOMO_FLOAT32 center[3];

  #if 0
  for (MESH_ELE_ID_TYPE v = 0; v < rp_printer_info.nVtx; v++)
    //size eup{
    *(Vtx + v * 3 + 0) = *(V0 + v * 3 + 0) * TOMO_FLOAT32(iVOXELFACTOR);
    *(Vtx + v * 3 + 1) = *(V0 + v * 3 + 1) * TOMO_FLOAT32(iVOXELFACTOR);
    *(Vtx + v * 3 + 2) = *(V0 + v * 3 + 2) * TOMO_FLOAT32(iVOXELFACTOR);
  }
  #endif

  //move center of mass to origin.
  AABB3Df.Set(printer_info.nVtx, V0);
  AABB3Df.GetCenter(center);
  for (TOMO_LONGINT  v = 0; v < printer_info.nVtx; v++)
  {
    *(Vtx + v * 3 + 0) = *(V0 + v * 3 + 0) - center[0];
    *(Vtx + v * 3 + 1) = *(V0 + v * 3 + 1) - center[1];
    *(Vtx + v * 3 + 2) = *(V0 + v * 3 + 2) - center[2];
  }

  //rotate
  SMatrix3f m3x3(printer_info.yaw, printer_info.pitch, printer_info.roll);
  for (TOMO_LONGINT v = 0; v < printer_info.nVtx; v++)
  {
    m3x3.Dot(&Vtx[v * 3], &Vtx[v * 3]);
    #ifdef _USE_VTX_NRM_FOR_PIXEL
    m3x3.Dot(&N0[v * 3], &Nrm[v * 3]);
    #endif
  }

 #ifndef _USE_VTX_NRM_FOR_PIXEL
  for (MESH_ELE_ID_TYPE t = 0; t < printer_info.nTri; t++)
 
  {
    m3x3.Dot(&N0[t * 3], &Nrm[t * 3]);
  }
#endif
 
  //translate so that corner lies on the origin.
  AABB3Df.Set(printer_info.nVtx, Vtx);
  for (TOMO_LONGINT v = 0; v < printer_info.nVtx; v++)
  {
    *(Vtx + v * 3 + 0) -= AABB3Df.x0;
    *(Vtx + v * 3 + 1) -= AABB3Df.y0;
    *(Vtx + v * 3 + 2) -= AABB3Df.z0;
  }

#ifdef _DEBUG
  AABB3Df.Set(printer_info.nVtx, Vtx);
#endif
}

void  STomoNV_Base::volToMass(void)
{
  vm_info.VolToMass(printer_info);
#if 0 
  vm_info.Vclad = printer_info.surface_area * printer_info.wall_thickness;
  vm_info.Vcore = vm_info.Vo - vm_info.Vclad;
  vm_info.Mclad = vm_info.Vclad * printer_info.Fclad * printer_info.PLA_density;
  vm_info.Mcore = vm_info.Vcore * printer_info.Fcore * printer_info.PLA_density;

#ifdef _EZAIR_THETAC_ZERO
  vm_info.Vss = vm_info.Vtc - vm_info.Vo;
#endif

  vm_info.Mo      = vm_info.Mcore + vm_info.Mclad;
  vm_info.Mss     = vm_info.Vss * printer_info.Fss* printer_info.Css * printer_info.PLA_density;
  vm_info.Mtotal  = vm_info.Mo + vm_info.Mss;
#endif
}
