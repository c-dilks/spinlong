#ifndef Environ_
#define Environ_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include "TSystem.h"

class Environ : public TObject
{
  public:
    Environ();

    Bool_t success;

    char SpinDir[512];
    char MassCutType[32];

    char MassCutsFile[1024];
    char ExclusionList[1024];

    Int_t year;
    char counts_file[512];
    char rtree_file[512];
    char pol_file[512];
    char output_dir[512];
    char redset_dir[512];

    Int_t which_eta_cut;

    Int_t PhiBins,EtaBins,PtBins,EnBins;
    Float_t PhiLow,PhiHigh;
    Float_t EtaLow,EtaHigh;
    Float_t PtLow,PtHigh;
    Float_t EnLow,EnHigh;

    Float_t PhiDiv(Int_t bin);
    Float_t EtaDiv(Int_t bin);
    Float_t EnDiv(Int_t bin);
    Float_t PtDiv(Int_t bin);
    Double_t Error(const char * bintype);

  protected:
    typedef std::map<Int_t,Float_t> DivMap;
    DivMap PhiDivMap,EtaDivMap,PtDivMap,EnDivMap;

    char MassCutsFileName[64];
    char ExclusionListName[64];

    ClassDef(Environ,1);
};

#endif
