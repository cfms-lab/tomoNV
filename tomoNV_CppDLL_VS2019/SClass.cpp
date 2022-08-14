#include "pch.h"
#include "SClass.h"

SClass::SClass()
{
  Init();
}

SClass::~SClass()
{
  Reset();
}

SClass::SClass(const SClass& Source)
{
  Init();
  _Copy(Source);
}

void	SClass::operator=(const SClass& Source)
{
  Reset();
  _Copy(Source);
}

void	SClass::_Copy(const SClass& Source)
{
  memcpy(iData, Source.iData, siData);
}

void	SClass::Reset(void)
{
  Init();
}

void	SClass::Init(void)
{
  memset( iData, 0x00, siData);
}

