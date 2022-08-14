#pragma once
#include "STPSlot.h"

using namespace Tomo;

class DLLEXPORT STomoNV_Base //base class of STomoNV_TMPxl class
{
public:
  STomoNV_Base();
  STomoNV_Base(const STomoNV_Base& Source);
  void	operator=(const STomoNV_Base& Source);
  void	_Copy(const STomoNV_Base& Source);
  ~STomoNV_Base();

  void	Reset(void);
  void	Init(void);

  S3DPrinterInfo    printer_info;
  STomoVolMassInfo  vm_info;
  TPSlotVector slotVec;
  STomoAABB2D AABB2D;

  void  Rotate(void); //rotate by (yaw, pitch, roll) and move to origin

  virtual void  Pixelize(void) {}
  virtual TPVector  slotsToPxls(enumPixelType _type) { TPVector tmp; return tmp;}//for rendering. time consuming. 
  virtual void      pxlsToSlots(TPVector& tri_pxls) {}
  virtual void  Pairing(void) {}//slot paring.
  virtual void  Calculate(void) {}//get Vss value
  void  volToMass(void);

  virtual TPVector GetSSPixels(void) { TPVector tmp; return tmp; }

};
