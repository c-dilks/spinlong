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
    TCUbits();

    void InitRPdefs();
    const char * RPname(Int_t idx0);
    Int_t RPidx(char * name);

    Bool_t Fired(char * trg);
    
    Bool_t FiredTOF();
    Bool_t FiredBBC();

    //Bool_t FiredRP(Int_t idx0);
    //Bool_t FiredRP(char * name0);

    Bool_t debug;
    Int_t NRP;

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

    std::map<std::string, Int_t> rp_idx; // rp trigger name to an index
    std::map<Int_t, std::string> rp_name; // rp index to trigger name

    ClassDef(TCUbits,1);
};

#endif
