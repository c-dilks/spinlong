// Implementation of L1 trigger bits (from TCU input)
// -- see ../TCUbits directory for further information

#ifndef TCUbits_
#define TCUbits_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>

#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TString.h"

#include "Environ.h"

class TCUbits : public TObject
{
  public:
    TCUbits(Environ * env0);

    Bool_t Fired(const char * trg);
    std::string WhichDSM(const char * trg);
    Int_t WhichTCUchan(const char * trg,std::string dsm_0);
    Int_t WhichBit(const char * trg);
    
    Bool_t debug;

    //-------------------------------
    // EVENT VARIABLES
    UInt_t lastdsm[8];
    //-------------------------------

  private:
    Environ * env;
    TTree * tcu_tr;
    TTree * tcuchan_tr;

    std::map<std::string, std::pair<std::string,UInt_t> > tcu_bit; // trigger --> (DSM,input bit)
    std::map<std::string, UInt_t> tcu_chan; // DSM --> TCU channel

    ClassDef(TCUbits,1);
};

#endif
