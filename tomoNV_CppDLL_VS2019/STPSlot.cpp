#include "pch.h"
#include "STPSlot.h"
#include <iostream> //cout. for debug
#include <math.h>//::pow()
#include "STomoPixel.h"//_weighted_addZ

using namespace Tomo;

void  splitAlBe_TMPxl(const TPVector& pxls, TPVector& al_pxls, TPVector& be_pxls)
{
  for (auto& pxl : pxls)
  {//we don't need pixels at the lateral side, i.e. nz==0.
    if      (pxl.nz * g_fNORMALFACTOR > g_fMARGIN * 10.)  {   al_pxls.push_back(pxl);   }
    else if (pxl.nz * g_fNORMALFACTOR < g_fMARGIN * -10.) {   be_pxls.push_back(pxl);   }
  }
}

TPVector STPSlot::createVoPxls(const TPVector& alpha, const TPVector& beta)
{
  TPVector Vo_pxls;//return data

  TPVector beta_minus;
  copyPxls(beta, beta_minus);
  negativeZ(beta_minus);

  TPVector temp_pxls;
  copyPxls(alpha, temp_pxls);
  temp_pxls.insert(temp_pxls.end(), beta_minus.begin(), beta_minus.end());
  if(temp_pxls.size() ==0) return Vo_pxls;//if none, return emtpy vector.

  STomoPixel Vo = _addZ(temp_pxls);
  Vo.z = _max(Vo.z , TOMO_INT16(0));//to reduce possible unexepcted noises.
  Vo.iTypeByte = typeVo;
  Vo_pxls.push_back(Vo);
  return Vo_pxls;
}

void _matchPairNumber_TMPxl(TPVector& al_slot, TPVector& be_slot)
{
  if (al_slot.size() > be_slot.size())
  {
    TPVector temp(al_slot.end() - be_slot.size(), al_slot.end()); al_slot = temp;//prefer Low Alpha.
  }
  else if (al_slot.size() < be_slot.size())
  {
    be_slot.resize(al_slot.size());//prefer High Beta. 
  }
}

TPVector createShadowAcceptor(const TPVector& nvB0, const TPVector& alpha0, bool _bUseExplicitSS, TPVector& explicit_ss_pxls)
{
  TPVector acceptor_pxls;
  for(auto& nvb: nvB0)
  {
#ifdef _DEBUG
TOMO_INT16 x = nvb.x;  TOMO_INT16 y = nvb.y;
#endif
    STomoPixel new_nbA = getNVA(nvb, alpha0);
    if(_bUseExplicitSS)
    {
      new_nbA.iTypeByte |= typeSSA;
      for (int z_offset = nvb.z; z_offset >= new_nbA.z; z_offset--)
      {
        STomoPixel new_ss(nvb);
        new_ss.z = z_offset;
        explicit_ss_pxls.push_back(new_ss);
      }
    }
    else
    {   new_nbA.iTypeByte |= typeNVA;    }
    acceptor_pxls.push_back(new_nbA);
  }
  return acceptor_pxls;
}

TPVector createVss_Explicit(const TPVector& ssb, const TPVector& ssa)
{
  TPVector vss_pxls;

  if(ssb.size() ==0) return vss_pxls;

  TPVector ssa_minus;
  copyPxls(ssa, ssa_minus);
  negativeZ(ssa_minus);

  TPVector pxls;
  pxls.insert(pxls.end(), ssb.begin(), ssb.end());
  pxls.insert(pxls.end(), ssa_minus.begin(), ssa_minus.end());

  if (pxls.size() == 0) return vss_pxls;//if none, return empty vector.

  STomoPixel Vss = _addZ(pxls);
  if (Vss.z != 0)//hide zero Vss pixels for pretty rendering.
  {
    Vss.iTypeByte = typeVss;
    vss_pxls.push_back(Vss);
  }

  return vss_pxls;
}
    

TPVector createVss_Implicit(
  const TPVector& alpha, const TPVector& beta,
  const TPVector& tc0,
  const TPVector& nvB, const TPVector& nvA0,
  std::function<TOMO_FLOAT32(TOMO_FLOAT32)> wFunc)
{
  TPVector vss_pxls;

  TPVector alpha_minus;
  copyPxls(alpha, alpha_minus);
  negativeZ(alpha_minus);

  TPVector nvB_minus;
  copyPxls(nvB, nvB_minus);
  negativeZ(nvB_minus);

  TPVector tc1;
  copyPxls(tc0, tc1);

  TPVector nva1;
  copyPxls(nvA0, nva1);

  zero_nZ(alpha_minus);// # we need only Nz of beta pixels
  zero_nZ(nvB_minus);
  zero_nZ(tc1);
  zero_nZ(nva1);

  TPVector pxls;
  pxls.insert(pxls.end(), alpha_minus.begin(), alpha_minus.end());
  pxls.insert(pxls.end(), beta.begin(), beta.end());
  pxls.insert(pxls.end(), tc1.begin(), tc1.end());
  pxls.insert(pxls.end(), nvB_minus.begin(), nvB_minus.end());
  pxls.insert(pxls.end(), nva1.begin(), nva1.end());

  if(pxls.size() ==0) return vss_pxls;//if none, return empty vector.

  STomoPixel Vss = _addZ(pxls, wFunc);
  if(Vss.z != 0)//hide zero Vss pixels for pretty rendering.
  {
    Vss.iTypeByte = typeVss;
    vss_pxls.push_back(Vss);
  }
  return vss_pxls;
}
//--------------------------------------------------------------------
STPSlot::STPSlot()
{
  Init();
}

STPSlot::~STPSlot()
{
  Reset();
}

STPSlot::STPSlot(const STPSlot& Source)
{
  Init();
  _Copy(Source);
}

void	STPSlot::operator=(const STPSlot& Source)
{
  Reset();
  _Copy(Source);
}

void	STPSlot::_Copy(const STPSlot& Source)
{
  pxls    = Source.pxls;
  vm_info = Source.vm_info;
  wFunc   = Source.wFunc;
  memcpy(iData, Source.iData, siData);
}

void	STPSlot::Reset(void)
{
  Init();
}

void	STPSlot::Init(void)
{
  pxls.clear();
  vm_info.Init();
  wFunc = nullptr;
  memset(iData, 0x00, siData);
}


void STPSlot::createAlBePxls(void)
{
  al_pxls.clear(); be_pxls.clear();
  sortSlotByZ(pxls, true);
  splitAlBe_TMPxl(pxls, al_pxls, be_pxls);
  removeZNearPxls(al_pxls);
  removeZNearPxls(be_pxls);
}

TPVector STPSlot::createTCPxls(void)
{
  if (al_pxls.size() > 0)
  {
    al_pxls.begin()->iTypeByte |= typeTC;
    TC_pxls.push_back(*al_pxls.begin());
  }
  return TC_pxls;
}

void  STPSlot::Pairing(TOMO_FLOAT32 theta_c, bool _bUseExplicitSS)
{
  if (pxls.size() == 0) return;

  //At the start, each slot has only alpha and beta pixels only.
  createAlBePxls();
  al_pxls.clear(); be_pxls.clear();
  sortSlotByZ(pxls, true);
  splitAlBe_TMPxl(pxls, al_pxls, be_pxls);
  removeZNearPxls(al_pxls);
  removeZNearPxls(be_pxls);

#ifdef _EZAIR_THETAC_ZERO
  createTCPxls();
  Vo_pxls = createVoPxls(al_pxls, be_pxls);
#else
  //Al-Be slot pairing.
  if (al_pxls.size() > 1) al_pxls = _matchAlternationWRT(al_pxls, be_pxls);
  if (be_pxls.size() > 1) be_pxls = _matchAlternationWRT(be_pxls, al_pxls);
  _matchPairNumber_TMPxl(al_pxls, be_pxls);

  Vo_pxls = createVoPxls(al_pxls, be_pxls);

  if (_bUseExplicitSS)//chapter 3.3.1 Explicit Method
  {
    SSB_pxls.clear(); SSA_pxls.clear(); SS_pxls.clear();
    SSB_pxls = createShadowCastorPixels(theta_c, be_pxls, true);
    SSB_pxls = removeZNearPxls(SSB_pxls);
    _matchPairWRT(SSB_pxls, be_pxls);

    SSA_pxls = createShadowAcceptor(SSB_pxls, al_pxls, true, SS_pxls);

    Vss_pxls = createVss_Explicit( SSB_pxls, SSA_pxls);//Eq. (17)
  }
  else//chapter 3.3.2 Implicit Method
  {
    createTCPxls();

    NVB_pxls.clear(); NVA_pxls.clear();
    NVB_pxls = createShadowCastorPixels(theta_c, be_pxls, false);//beta pixles have been changed in matchAlBeSlots(). remake them.
    NVB_pxls = removeZNearPxls(NVB_pxls);
    _matchPairWRT(NVB_pxls, be_pxls);

    NVA_pxls = createShadowAcceptor(NVB_pxls, al_pxls, false, SS_pxls);
     
    Vss_pxls = createVss_Implicit(al_pxls, be_pxls, TC_pxls, NVB_pxls, NVA_pxls, wFunc);//Eq. (19)
  }
#endif
  pxls.clear();
  pxls << al_pxls << be_pxls;
  pxls << TC_pxls << NVB_pxls << NVA_pxls;
  pxls << SSB_pxls << SSA_pxls << SS_pxls;
  pxls << Vo_pxls << Vss_pxls;
  sortSlotByZ(pxls, true);
}

void  STPSlot::Calculate(void)//calculate local vol. & mass.
{
  vm_info.Init();
  for (auto& pxl : al_pxls)    {     vm_info.Va  += pxl.z;    }
  for (auto& pxl : be_pxls)    {     vm_info.Vb  += pxl.z;    }
  for (auto& pxl : TC_pxls)    {     vm_info.Vtc += pxl.z;    }
  for (auto& pxl : NVB_pxls)   {     vm_info.Vnv += pxl.z;    }
  for (auto& pxl : NVA_pxls)   {     vm_info.Vnv -= pxl.z;    }
#ifdef _DEBUG
  for (auto& pxl : SSB_pxls)   {     vm_info.SS_vol += pxl.z; }
  for (auto& pxl : SSA_pxls)   {     vm_info.SS_vol -= pxl.z; }
#endif
        
  if (Vo_pxls.size() > 0)
  {
    vm_info.Vo = TOMO_FLOAT32(Vo_pxls.begin()->z); //Eq. (4)
  }
  if (Vss_pxls.size() > 0)
  {
    vm_info.Vss = TOMO_FLOAT32(Vss_pxls.begin()->z); // Eq. (5)
  }


}

TPVector  STPSlot::slotToPxls(enumPixelType _type)
{
  TPVector pxls_with_types;

  int iTypeBit = toTypeByte(_type);

  for (auto& pxl : pxls)
  {
    if (pxl.iTypeByte & iTypeBit)
    {
      pxls_with_types.push_back(pxl);
    }
  }

  return pxls_with_types;
}

#if 1
TPVector  STPSlot::GetSSPxls(bool _bUseExplicitSS)//prepare support structure pxls for rendering.
{
  if (_bUseExplicitSS)
  {
    return SS_pxls;
  }

  if (Vss_pxls.size() == 0 || Vss_pxls.begin()->z <= 0) return SS_pxls;

  SS_pxls.clear();
  int ss_start_z = -1, ss_end_z = -1;
  size_t slot_len = pxls.size();
  for (TPIterator pIt = pxls.begin() ; pIt != pxls.end() ; pIt++)
  {
    if ((pIt->iTypeByte & typeBe) && !(pIt->iTypeByte & typeNVB))
    {//start of SS pxls
      ss_start_z = pIt->z;
    }

    if ((ss_start_z>=0) && (pIt->iTypeByte & typeAl))
    {
      ss_end_z = pIt->z;
    }

    if ((ss_start_z >= 0) && (pIt == pxls.end()- 1/*right before bottom plate*/))
    {
      ss_end_z = 0;
    }

    if (ss_start_z > -1 && ss_end_z > -1)
    {
      for (int ss_z = ss_start_z; ss_z >= ss_end_z; ss_z--)
      {
        STomoPixel new_pxl(pIt->x, pIt->y, ss_z, 0, 0, 0);
        SS_pxls.push_back(new_pxl);
      }
      ss_start_z = -1; ss_end_z = -1;//reset
    }
  }

  return SS_pxls;
}
#endif
