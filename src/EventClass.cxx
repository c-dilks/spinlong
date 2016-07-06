#include "EventClass.h"

ClassImp(EventClass)

namespace
{
  const Int_t N_CLASS = 5;
  enum class_enum
  {
    kSph,
    kPi0,
    kThr,
    kEtm,
    kDpi
  };

  const Float_t Zd = 720;  // distance from IP [cm]
  const Float_t Sd = 3.8; // small cell dimension [cm]
  const Float_t Ld = 5.8; // large cell dimension [cm]

  const Float_t pi0_mass = 0.135; // pi0 mass [GeV]
  const Float_t etm_mass = 0.548; // eta mass [GeV]
  //const Float_t jps_mass = 3.097; // j/psi mass [GeV]
};

// constructor
// DoNotInitKinBounds will not initialise KinBounds; this is only used in 
// add_diag.C, because add_diag.C produces threshtr, a dependency of KinBounds
EventClass::EventClass(Environ * env0, Bool_t DoNotInitKinBounds)
{
  N = N_CLASS;

  // event class names
  class_name.insert(std::pair<Int_t,char*>(kSph,"sph")); // single photons
  class_name.insert(std::pair<Int_t,char*>(kPi0,"pi0")); // pi0's
  class_name.insert(std::pair<Int_t,char*>(kThr,"thr")); // three or more photons
  class_name.insert(std::pair<Int_t,char*>(kEtm,"etm")); // eta's
  class_name.insert(std::pair<Int_t,char*>(kDpi,"dpi")); // di-pi0's

  // event class titles
  class_title.insert(std::pair<Int_t,char*>(kSph,"Single #gamma")); // single photons
  class_title.insert(std::pair<Int_t,char*>(kPi0,"#pi^{0}")); // pi0's
  class_title.insert(std::pair<Int_t,char*>(kThr,"N_{#gamma}>2")); // three or more photons
  class_title.insert(std::pair<Int_t,char*>(kEtm,"#eta-meson")); // eta's
  class_title.insert(std::pair<Int_t,char*>(kDpi,"di-#pi{0}")); // di-pi0's

  // event class index from name
  for(Int_t n=0; n<N_CLASS; n++)
    class_idx.insert(std::pair<std::string,Int_t>(std::string(class_name[n]),n));

  // read kinematic-dependent mass cuts
  env = env0;
  mass_tr = new TTree();
  char mass_tr_file[512];
  sprintf(mass_tr_file,"%s/mass_cuts.dat",env->SpinDir);
  mass_tr->ReadFile(env->MassCutsFile,"kbinL/F:kbinH/F:massL/F:massM/F:massH/F");
  mass_tr->SetBranchAddress("kbinL",&kbinL);
  mass_tr->SetBranchAddress("kbinH",&kbinH);
  mass_tr->SetBranchAddress("massL",&massL);
  mass_tr->SetBranchAddress("massM",&massM);
  mass_tr->SetBranchAddress("massH",&massH);

  // import exclusion list
  exclude_tr = new TTree();
  char exclude_tr_file[512];
  sprintf(exclude_tr_file,"%s/exclusion_list",env->SpinDir);
  exclude_tr->ReadFile(env->ExclusionList,"exc_run/I");
  exclude_tr->SetBranchAddress("exc_run",&exc_run);

  runnum_tmp=0;
  exclude_run=false;

  if(!DoNotInitKinBounds) 
    KB = new KinBounds(env);
};


// sets kinematic variables for the event, used to check cuts with Valid()
void EventClass::SetKinematics(Int_t runnum_,
                               Float_t E12_,
                               Float_t Pt_,
                               Float_t Eta_,
                               Float_t Phi_,
                               Float_t M12_,
                               Float_t Z_,
                               Float_t N12_,
                               Int_t ClIndex_)
{
  runnum = runnum_;
  E12 = E12_;
  Pt = Pt_;
  Eta = Eta_;
  Phi = Phi_;
  M12 = M12_;
  Z = Z_;
  N12 = N12_;
  ClIndex = ClIndex_;
};

// appends kinematics to "kinematics" list: kinematics for each cluster;
// used for multi-hit events which aren't simply two photons
void EventClass::AppendKinematics(Float_t E12_,
                                  Float_t Pt_,
                                  Float_t Eta_,
                                  Float_t Phi_,
                                  Float_t M12_,
                                  Float_t Z_,
                                  Float_t N12_,
                                  Int_t ClIndex_)
{
  E12list[ClIndex_] = E12_;
  Ptlist[ClIndex_] = Pt_;
  Etalist[ClIndex_] = Eta_;
  Philist[ClIndex_] = Phi_;
  M12list[ClIndex_] = M12_;
  Zlist[ClIndex_] = Z_;
  N12list[ClIndex_] = N12_;
};


// copies current kinematics under use to "tmp" kinematics;
// RecallKinematics does the opposite
void EventClass::StoreKinematics()
{
  E12tmp = E12;
  Pttmp = Pt;
  Etatmp = Eta;
  Phitmp = Phi;
  M12tmp = M12;
  Ztmp = Z;
  N12tmp = N12;
  ClIndextmp = ClIndex;
};
void EventClass::RecallKinematics()
{
  E12 = E12tmp;
  Pt = Pttmp;
  Eta = Etatmp;
  Phi = Phitmp;
  M12 = M12tmp;
  Z = Ztmp;
  N12 = N12tmp;
  ClIndex = ClIndextmp;
};


// returns event class index
Int_t EventClass::Idx(char * name)
{
  Int_t retval;
  try { retval = class_idx.at(std::string(name)); }
  catch(const std::out_of_range& e) {
    retval=-1; 
    //fprintf(stderr,"EventClass::Idx out of range\n");
  };
  return retval;
};

// returns event class name, given index
char * EventClass::Name(Int_t idx)
{
  if(idx>=0 && idx<N) return class_name.at(idx);
  else return "";
};

// returns event class title, given index
char * EventClass::Title(Int_t idx)
{
  if(idx>=0 && idx<N) return class_title.at(idx);
  else return "";
};

// returns event class title, given event class name
char * EventClass::Title(char * name)
{
  return class_title.at(Idx(name));
};


//////////////////////////
//   EVENT CLASS CUTS   //
//////////////////////////
// returns true if the event passes the cuts
// -- if trig_index>=0, we check 
//    tighter cuts on kinematics via KinBounds; this is useful
//    for selecting events for A_LL
// -- if trig_index<0, we just check kinematics are in ranges specified
//    in Bin_Splitter.C; this is
//    useful for KinVarVsRun.C and DiagnosticsOne.C, where we'd like
//    to see the full distribution
// -- see doc_diagram.pdf for a diagram of the flow of information through
//    all of these classes
//
Bool_t EventClass::Valid(Int_t idx, Int_t trig_index)
{
  // check if on exclusion list
  if(runnum!=runnum_tmp)
  {
    runnum_tmp = runnum;
    exclude_run = ExcludedRun();

    /*
    exclude_run=false;
    for(Int_t e=0; e<exclude_tr->GetEntries(); e++)
    {
      exclude_tr->GetEntry(e);
      if(runnum==exc_run) exclude_run=true;
    };
    */
  };
  if(exclude_run) return false;


  Bool_t validity = false;


  // single photon cuts
  if(idx==kSph)
  {
    if( fabs(N12-1)<0.1 /*&&*/
        /*E12>10 &&*/
        /*Pt>1*/ ) validity = true;
  }

  // pi0 cuts
  else if(idx==kPi0)
  {
    if( fabs(N12-2)<0.1 &&
        ClIndex==0 &&
        Z<0.8 &&
        CheckMass(M12) /*&&*/
        /*E12>10 &&*/
        /*Pt>1*/ ) validity = true;
  }

  // three or more photons cuts
  else if(idx==kThr)
  {
    if( N12>2 &&
        M12>0.7 /*&& */
        /*E12>10 &&*/
        /*Pt>1*/ ) validity = true;
  }

  // eta meson cuts
  else if(idx==kEtm)
  {
    if( fabs(N12-2)<0.1 &&
        Z<0.8 &&
        fabs(M12-etm_mass)<0.15 &&
        FiducialGeom(Eta,Phi,1.5) /*&&*/
        /*E12>10 &&*/
        /*Pt>1*/ ) validity = true;
  }

  // j/psi cuts
  /*
  else if(idx==kJps)
  {
    if( fabs(N12-2)<0.1 &&
        Z<0.3 &&
        M12>2 && M12<3.5 &&
        Eta>3.2 &&
        FiducialGeom(Eta,Phi,1.5) &&
        E12>60 && E12<100 &&
        Pt>1) validity = true;
  }
  */

  else if(idx==kDpi) {
    StoreKinematics();
    picnt=0;
    //printf("looping over %d clusters:\n",Nclust);
    for(int clu=0; clu<Nclust; clu++) {
      if(picnt<2) {
        // set main kinematic variables to those of this cluster
        SetKinematics(runnum,
                      E12list[clu],
                      Ptlist[clu],
                      Etalist[clu],
                      Philist[clu],
                      M12list[clu],
                      Zlist[clu],
                      N12list[clu],
                      0);
        // see if this cluster contains a pi0
        //if(Valid(kPi0,trig_index)) {
        // see if this first cluster contains pi0 and second cluster is anything E>10 GeV
        if(Valid(kPi0,trig_index) || (picnt==1 && E12list[clu]>10.0)) {
          dipi_E12[picnt] = E12;
          dipi_Pt[picnt] = Pt;
          dipi_Eta[picnt] = Eta;
          dipi_Phi[picnt] = Phi;
          dipi_M12[picnt] = M12;
          dipi_Z[picnt] = Z;
          picnt++;
          //printf("%d ---%s\n",clu,(picnt==1)?"pion":"second");
        };
      };
    };
    if(picnt==2) {
      // compute Z-component of cross product of momenta
      delta_phi = dipi_Phi[1]-dipi_Phi[2];
      sin_delta_phi = sin(delta_phi);
      pt_prod = dipi_Pt[0]*dipi_Pt[1];
      cross_prod_z = pt_prod * sin_delta_phi;
      validity = true;
    };
    RecallKinematics();
  };
   
   // ------ //


  // check tighter kinematic cuts (run-by-run, trigger-by-trigger)
  if(validity && idx!=kDpi) {
    if(trig_index>=0) {
      validity=false;
      if(KB->PtInRange(Pt,runnum,idx,trig_index) 
         /*&& KB->EnInRange(E12,runnum,idx,trig_index)*/) {
        validity = true;
      };
    };
  };


  return validity;
};
//////////////////////////


// checks Valid(), but ignores Mass cut
Bool_t EventClass::ValidWithoutMcut(Int_t idx)
{
  Float_t M12_tmp = M12; // store M12 to tmp variable
  Bool_t boole;
  // set M12 to optimal meson masses
  if(idx==kPi0) M12 = pi0_mass;
  else if(idx==kEtm) M12 = etm_mass;
  /*else if(idx==kJps) M12 = jps_mass;*/
  boole = Valid(idx);
  M12 = M12_tmp; // restore value of M12
  return boole;
};


// checks Valid(), but ignores Z cut
Bool_t EventClass::ValidWithoutZcut(Int_t idx)
{
  Float_t Z_tmp = Z; // store Z to tmp variable
  Float_t boole;
  Z = 0; // set Z to optimal value
  boole = Valid(idx);
  Z = Z_tmp; // restore value of Z
  return boole;
};


// check E-dependent or Pt-dependent mass cut; returns true if passed
Bool_t EventClass::CheckMass(Float_t M12_)
{
  for(Int_t q=0; q<mass_tr->GetEntries(); q++)
  {
    mass_tr->GetEntry(q);
    if( (M12>=massL && M12<=massH) &&
        ( (!strcmp(env->MassCutType,"en") && E12>=kbinL && E12<=kbinH) ||
          (!strcmp(env->MassCutType,"pt") && Pt>=kbinL && Pt<=kbinH) )) return true;
  };
  return false;
};


// checks if Eta and Phi are within fiducial area cut, defined as Cd-cells
// from cell boundaries
Bool_t EventClass::FiducialGeom(Float_t Eta_, Float_t Phi_, Float_t Cd)
{
  Bool_t boole = false;
  Theta = 2*atan2(exp(-1*Eta_),1); // polar angle

  Xd =  Zd * tan(Theta) * cos(Phi_);
  Yd =  Zd * tan(Theta) * sin(Phi_);

  // small cell fiducial area
  if( ( fabs(Xd)<(12-Cd)*Sd && fabs(Yd)<(12-Cd)*Sd ) && 
     !( fabs(Xd)<( 5+Cd)*Sd && fabs(Yd)<( 5+Cd)*Sd )) boole = true;

  // large cell fiducial area
  if( ( fabs(Xd)<(17-Cd)*Ld && fabs(Yd)<(17-Cd)*Ld ) &&
     !( fabs(Xd)<( 8+Cd)*Ld && fabs(Yd)<( 8+Cd)*Ld ) &&
     fabs(Xd+Yd)<(26-sqrt(2)*Cd)*Ld &&
     fabs(Xd-Yd)<(26-sqrt(2)*Cd)*Ld) boole = true;


  // uncomment these lines for cut which was used for spin2014...
  //if(Cd!=0) boole = false;
  //if(Cd!=0 && Eta_>=2.5 && Eta_<=4) boole = true;

  return boole;
};


// returns true if run is on exclusion list
Bool_t EventClass::ExcludedRun() {
  for(Int_t e=0; e<exclude_tr->GetEntries(); e++)
  {
    exclude_tr->GetEntry(e);
    if(runnum==exc_run) return true;
  };
  return false;
};
