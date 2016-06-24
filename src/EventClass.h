#ifndef EventClass_
#define EventClass_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"

#include "Environ.h"
#include "KinBounds.h"


class EventClass : public TObject
{
  public:
    EventClass(Environ * env0, Bool_t DoNotInitKinBounds=false);
    void SetKinematics(Int_t runnum_,
                       Float_t E12_,
                       Float_t Pt_,
                       Float_t Eta_,
                       Float_t Phi_,
                       Float_t M12_,
                       Float_t Z_,
                       Float_t N12_,
                       Int_t ClIndex_);
    Int_t Idx(char * name);
    char * Name(Int_t idx);
    char * Title(Int_t idx);
    char * Title(char * name);
    Bool_t Valid(Int_t idx,Int_t trig_index=-1); // returns true if cuts pass
    Bool_t ValidWithoutMcut(Int_t idx); // Valid(), but don't cut on mass
    Bool_t ValidWithoutZcut(Int_t idx); // Valid(), but don't cut on Z
    Bool_t CheckMass(Float_t M12_);
    Bool_t FiducialGeom(Float_t Eta_,Float_t Phi_, Float_t Cd);
    Bool_t ExcludedRun(); // true if excluded

    Int_t N;
    Int_t runnum,ClIndex;
    Float_t E12,Pt,Eta,Phi,M12,Z,N12;

    Float_t Theta,Xd,Yd;

  protected:
    std::map<Int_t, char*> class_name; // idx --> name
    std::map<std::string, Int_t> class_idx;  // name --> idx
    std::map<Int_t, char*> class_title; // idx --> title

    Int_t runnum_tmp;
    Bool_t exclude_run;

    Environ * env;
    KinBounds * KB;
    TTree * mass_tr;
    Float_t kbinL,kbinH,massL,massM,massH;
    TTree * exclude_tr;
    Int_t exc_run;

  ClassDef(EventClass,1);
};

#endif
