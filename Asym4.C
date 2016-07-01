// computes asymmetry between fully summed (i.e. over all runs) phi distributions
//    (see staszak thesis, eq. 7.5)
//
//  - this script does all asymmetries; currently A_LL and A_L for both beams, each classified
//    by an "asymmetry number", defined as the Rellum number used to compute the asymmetry:
//    -- asym=1 :: A_L yellow
//    -- asym=2 :: A_L blue
//    -- asym=3 :: A_LL double helicity asymmetry
//  
//  - filter type: defines a filter for the data, useful for consistency checks
//    -- all: entire data set
//    -- run: only runs from filter_low to filter_high
//      -- runout: only runs excluding those from filter_low to filter_high
//      -- runeven: only even run numbers
//      -- runodd: only odd run numbers
//    -- fill: only fills from filter_low to filter_high
//  --> use asym_call* scripts to call different filters
//

void Asym4(const char * evclass="pi0", const char * filter_type="all",Int_t filter_low=0, Int_t filter_high=0)
{
  // ANALYSIS TYPE -- use this swtich to change between a longitudinal and a transverse analysis
  enum atypes {kLong,kTrans};
  const Int_t ANALYSIS_TYPE = kLong;  // ---- SWITCH ---- //

  Int_t NPARAM_tmp;
  switch(ANALYSIS_TYPE) {
    case kLong:
      NPARAM_tmp = 1;
      break;
    case kTrans:
      NPRAM_tmp = 2;
      break;
  };
  const Int_t NPARAM = NPARAM_tmp;


  const Float_t pi=3.1415;
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  TString infile_n = Form("%s/all.root",RD->env->phiset_dir);

  TFile * infile = new TFile(infile_n.Data(),"READ");

  // get bins from environment
  Int_t phi_bins0 = RD->env->PhiBins; const Int_t phi_bins = phi_bins0;
  Int_t eta_bins0 = RD->env->EtaBins; const Int_t eta_bins = eta_bins0;
  Int_t en_bins0 = RD->env->EnBins; const Int_t en_bins = en_bins0;
  Int_t pt_bins0 = RD->env->PtBins; const Int_t pt_bins = pt_bins0;



  // read TObjArrays of phi distributions
  TObjArray * phi_dist_arr[4][eta_bins][pt_bins][en_bins];
  char phi_dist_arr_n[4][eta_bins][pt_bins][en_bins][64];
  Int_t NRUNS_tmp=0;
  Int_t ARR_size;
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(phi_dist_arr_n[s][g][p][e],"%s/phi_dist_%s_s%d_g%d_p%d_e%d",evclass,evclass,s,g,p,e);
          phi_dist_arr[s][g][p][e] = (TObjArray*) infile->Get(phi_dist_arr_n[s][g][p][e]);
          printf("phi_dist_arr[%d][%d][%d][%d] @ %p\n",s,g,p,e,(void*)phi_dist_arr[s][g][p][e]);
          if(s==0 && g==0 && p==0 && e==0)
          {
            ARR_size=phi_dist_arr[s][g][p][e]->GetEntries();
            for(Int_t kk=0; kk<ARR_size; kk++)
            {
              if(((TH1D*)(phi_dist_arr[s][g][p][e]->At(kk)))->GetEntries() > 0) NRUNS_tmp++;
            }
          }
          else
          {
            if(phi_dist_arr[s][g][p][e]->GetEntries() != ARR_size)
            {
              fprintf(stderr,"ERROR: TObjArrays have different sizes\n");
              return;
            };
          };
        };
      };
    };
  };
  const Int_t NRUNS = NRUNS_tmp;
  printf("ARR_size=%d\n",ARR_size);
  printf("NRUNS=%d\n",NRUNS);



  // build summed phi distributions with appropriate polarization and rellum weights
  // *_num = numerator in MLM calculation
  // *_den = denominator in MLM calculation
  // *_num_e = numerator in statistical error
  // *_den_e = denominator in statistical error
  const Int_t asym_bins=4;
  TH1D * phi_dist_num[asym_bins][4][eta_bins][pt_bins][en_bins]; // [asymmetry number] [spinbit] [eta] [pt] [en] 
  TH1D * phi_dist_den[asym_bins][4][eta_bins][pt_bins][en_bins];
  TH1D * phi_dist_num_e[asym_bins][4][eta_bins][pt_bins][en_bins]; 
  TH1D * phi_dist_den_e[asym_bins][4][eta_bins][pt_bins][en_bins];
  char phi_dist_num_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_den_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_num_e_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  char phi_dist_den_e_n[asym_bins][4][eta_bins][pt_bins][en_bins][128];
  Double_t yield[4][eta_bins][pt_bins][en_bins]; // counts yields (only for a==3); used for error analysis
  Int_t runnum;
  Float_t rellum,polar_b,polar_y;
  Float_t weight_num,weight_den; // weight for asymmetry (MLM for A_LL)
  Float_t weight_num_e,weight_den_e; // weight for statistical error (by-product of MLM for A_LL);
  Int_t fill,pattern;
  Bool_t isConsistent;
  TString tmpstr;
  gROOT->ProcessLine(".! touch runlist; rm runlist; touch runlist");
  for(Int_t s=0; s<4; s++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          // initialise arrays
          yield[s][g][p][e] = 0;
        };
      };
    };
  };
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            if(a==3) yield[s][g][p][e] = 0; // initialise yield counter

            sprintf(phi_dist_num_n[a][s][g][p][e],"phi_num_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            sprintf(phi_dist_den_n[a][s][g][p][e],"phi_den_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            sprintf(phi_dist_num_e_n[a][s][g][p][e],"phi_num_e_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            sprintf(phi_dist_den_e_n[a][s][g][p][e],"phi_den_e_a%d_s%d_g%d_p%d_e%d",a,s,g,p,e);
            phi_dist_num[a][s][g][p][e] = new TH1D(phi_dist_num_n[a][s][g][p][e],phi_dist_num_n[a][s][g][p][e],
              phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
            phi_dist_den[a][s][g][p][e] = new TH1D(phi_dist_den_n[a][s][g][p][e],phi_dist_den_n[a][s][g][p][e],
              phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
            phi_dist_num_e[a][s][g][p][e] = new TH1D(phi_dist_num_e_n[a][s][g][p][e],phi_dist_num_e_n[a][s][g][p][e],
              phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
            phi_dist_den_e[a][s][g][p][e] = new TH1D(phi_dist_den_e_n[a][s][g][p][e],phi_dist_den_e_n[a][s][g][p][e],
              phi_bins,RD->env->PhiLow,RD->env->PhiHigh);

            for(Int_t r=0; r<ARR_size; r++)
            {
              tmpstr = TString(phi_dist_arr[0][g][p][e]->At(r)->GetName());
              sscanf(tmpstr(tmpstr.Length()-8,tmpstr.Length()).Data(),"%d",&runnum);

              //printf("phi_dist_name=%s\n",phi_dist_arr[0][g][p][e]->At(r)->GetName());
              //printf("evclass=%s runnum=%d\n",evclass,runnum);

              // RELLUM SELECTION HERE
              rellum = RD->Rellum(runnum,a,"vpd"); // note that asym no. = rellum no. needed for this asymmetry
              //rellum=1;     // for testing
              polar_b = RD->BluePol(runnum);
              polar_y = RD->YellPol(runnum);
              fill = RD->GetFill(runnum);
              isConsistent = RD->RellumConsistent(runnum);
              pattern = RD->Pattern(runnum);


              // polarization and rellum weighting; depends on which asymmetry
              if(a==3)
              {
                weight_num = polar_b * polar_y;
                weight_den = pow(polar_b * polar_y, 2);
                weight_num_e = pow(polar_b * polar_y, 2);
                weight_den_e = pow(polar_b * polar_y, 2);
                if(s==1 || s==2) 
                {
                  weight_num *= rellum;
                  weight_den *= rellum;
                  weight_num_e *= pow(rellum, 2);
                  weight_den_e *= rellum;
                };
              }
              else if(a==1)
              {
                weight_num = polar_y;
                weight_den = pow(polar_y, 2);
                weight_num_e = pow(polar_y, 2);
                weight_den_e = pow(polar_y, 2);
                if(s==0 || s==2) 
                {
                  weight_num *= rellum;
                  weight_den *= rellum;
                  weight_num_e *= pow(rellum, 2);
                  weight_den_e *= rellum;
                };
              }
              else if(a==2)
              {
                weight_num = polar_b;
                weight_den = pow(polar_b, 2);
                weight_num_e = pow(polar_b, 2);
                weight_den_e = pow(polar_b, 2);
                if(s==0 || s==1) 
                {
                  weight_num *= rellum;
                  weight_den *= rellum;
                  weight_num_e *= pow(rellum, 2);
                  weight_den_e *= rellum;
                };
              };

              
              // print out runlist with fill no. and R3
              if(a==3 && s==0 && g==0 && p==0 && e==0 && ((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)))->GetEntries()>0)
              {
                gSystem->RedirectOutput("runlist");
                printf("%d %d %d %f\n",r,runnum,fill,rellum);
                gSystem->RedirectOutput(0);
              };

              // add weighted dists if we have consistent rellum measurement and filter_type 
              // defined above is passed
              if(isConsistent && ((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)))->GetEntries() > 0)
              {
                if( ( !strcmp(filter_type,"fill") && (fill>=filter_low && fill<=filter_high) ) ||
                    ( !strcmp(filter_type,"run") && (runnum>=filter_low && runnum<=filter_high) ) ||
                    ( !strcmp(filter_type,"runout") && !(runnum>=filter_low && runnum<=filter_high) ) ||
                    ( !strcmp(filter_type,"runeven") && (runnum % 2)==0 ) ||
                    ( !strcmp(filter_type,"runodd") && (runnum % 2)==1 ) ||
                      !strcmp(filter_type,"all"))
                {
                  phi_dist_num[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_num);
                  phi_dist_den[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_den);
                  phi_dist_num_e[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_num_e);
                  phi_dist_den_e[a][s][g][p][e]->Add((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)),weight_den_e);
                  if(a==3) yield[s][g][p][e] += ((TH1D*)(phi_dist_arr[s][g][p][e]->At(r)))->GetEntries(); // increment yield counter
                };
              };
            };
          };
        };
      };
    };
  };


  // compute asymmetry; polarization and rellum have been corrected for above
  // asym = (ll_num - rr_num) / (ll_den + rr_den)
  // -- a=1 :: ll=s1+s3 :: rr=s0+s2
  // -- a=2 :: ll=s2+s3 :: rr=s0+s1
  // -- a=3 :: ll=s0+s3 :: rr=s1+s2
  // also compute statistical error bars
  // asym_e = sqrt(ll_num_e + rr_num_e) / (ll_den_e + rr_den_e)
  TH1D * dist_ll_num[asym_bins][eta_bins][pt_bins][en_bins]; // left & right asym terms
  TH1D * dist_rr_num[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_ll_den[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_rr_den[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_ll_num_e[asym_bins][eta_bins][pt_bins][en_bins]; // left & right error terms
  TH1D * dist_rr_num_e[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_ll_den_e[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * dist_rr_den_e[asym_bins][eta_bins][pt_bins][en_bins]; 
  TH1D * numer[asym_bins][eta_bins][pt_bins][en_bins]; // numer = ll_num - rr_num
  TH1D * denom[asym_bins][eta_bins][pt_bins][en_bins]; // denom = ll_den + rr_den
  TH1D * numer_e[asym_bins][eta_bins][pt_bins][en_bins]; // numer_e = ll_num_e + rr_num_e
  TH1D * denom_e[asym_bins][eta_bins][pt_bins][en_bins]; // denom_e = ll_den_e + rr_den_e
  TH1D * numer_e_sqrt[asym_bins][eta_bins][pt_bins][en_bins]; // numer_e_sqrt = sqrt(numer_e)
  TH1D * asym[asym_bins][eta_bins][pt_bins][en_bins]; // asym = numer / denom
  TH1D * asym_e[asym_bins][eta_bins][pt_bins][en_bins]; // asym_e = numer_e_sqrt / denom_e
  Int_t asym_pts[asym_bins][eta_bins][pt_bins][en_bins]; // number of points in asym

  char dist_ll_num_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_num_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_ll_den_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_den_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_ll_num_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_num_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_ll_den_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char dist_rr_den_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char numer_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char denom_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char numer_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char denom_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char numer_e_sqrt_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char asym_n[asym_bins][eta_bins][pt_bins][en_bins][128];
  char asym_t[asym_bins][eta_bins][pt_bins][en_bins][256];
  char asym_e_n[asym_bins][eta_bins][pt_bins][en_bins][128];

  Float_t p0,p0e,chi2,ndf;
  Float_t bc[4];
  Float_t bcent;
  Double_t bc_e,bc_ee;
  Int_t runnum_0;
  Int_t qq;
  Float_t asym_value[asym_bins][eta_bins][pt_bins][en_bins];
  Int_t asym_value_cnt[asym_bins][eta_bins][pt_bins][en_bins];
  TF1 * asym_fit[asym_bins][eta_bins][pt_bins][en_bins];
  char asym_fit_n[asym_bins][eta_bins][pt_bins][en_bins][32];
  Float_t asym_tmp;
  Float_t asym_max[asym_bins][eta_bins][pt_bins][en_bins];
  Float_t asym_min[asym_bins][eta_bins][pt_bins][en_bins];
  char var_str[16]; strcpy(var_str,"#phi");


  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          asym_value[a][g][p][e]=0.0;
          asym_value_cnt[a][g][p][e]=0;
          asym_max[a][g][p][e]=0;
          asym_min[a][g][p][e]=0;
        };
      };
    };
  };

  char asym_title[asym_bins][60];
  if(ANALYSIS_TYPE==kTrans) {
    strcpy(asym_title[1],"Transverse SSA (Y)");
    strcpy(asym_title[2],"Transverse SSA (B)");
    strcpy(asym_title[3],"Transverse DSA");
  }
  else if(ANALYSIS_TYPE==kLong) {
    strcpy(asym_title[1],"Longitudinal SSA (Y)");
    strcpy(asym_title[2],"Longitudinal SSA (B)");
    strcpy(asym_title[3],"Longitudinal DSA");
  };

  char asym_title_kd[NPARAM][asym_bins][20]; // for kin dep plots
  if(ANALYSIS_TYPE==kTrans) {
    strcpy(asym_title_kd[0][1],"R (Y)");
    strcpy(asym_title_kd[0][2],"R (B)");
    strcpy(asym_title_kd[0][3],"A_{#Sigma}");
    strcpy(asym_title_kd[1][1],"A_{N}^{Y}");
    strcpy(asym_title_kd[1][2],"A_{N}^{B}");
    strcpy(asym_title_kd[1][3],"A_{TT}");
  }
  else if(ANALYSIS_TYPE==kLong) {
    strcpy(asym_title_kd[0][1],"A_{L}^{Y}");
    strcpy(asym_title_kd[0][2],"A_{L}^{B}");
    strcpy(asym_title_kd[0][3],"A_{LL}");
  };

  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          // initialise ll,rr,numer,denom,asym histograms
          sprintf(dist_ll_num_n[a][g][p][e],"dist_ll_num_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_num_n[a][g][p][e],"dist_rr_num_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_ll_den_n[a][g][p][e],"dist_ll_den_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_den_n[a][g][p][e],"dist_rr_den_a%d_g%d_p%d_e%d",a,g,p,e);

          sprintf(dist_ll_num_e_n[a][g][p][e],"dist_ll_num_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_num_e_n[a][g][p][e],"dist_rr_num_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_ll_den_e_n[a][g][p][e],"dist_ll_den_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(dist_rr_den_e_n[a][g][p][e],"dist_rr_den_e_a%d_g%d_p%d_e%d",a,g,p,e);

          sprintf(numer_n[a][g][p][e],"numer_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(denom_n[a][g][p][e],"denom_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(numer_e_n[a][g][p][e],"numer_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(denom_e_n[a][g][p][e],"denom_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(numer_e_sqrt_n[a][g][p][e],"numer_e_sqrt_a%d_g%d_p%d_e%d",a,g,p,e);

          sprintf(asym_n[a][g][p][e],"asym_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(asym_e_n[a][g][p][e],"asym_e_a%d_g%d_p%d_e%d",a,g,p,e);
          sprintf(asym_t[a][g][p][e],
             "%s vs. %s :: #eta#in[%.2f,%.2f), p_{T}#in[%.2f,%.2f), E#in[%.2f,%.2f) (runsum)",
             asym_title[a],var_str,
             RD->env->EtaDiv(g),RD->env->EtaDiv(g+1),
             RD->env->PtDiv(p),RD->env->PtDiv(p+1),
             RD->env->EnDiv(e),RD->env->EnDiv(e+1));


          dist_ll_num[a][g][p][e] = new TH1D(dist_ll_num_n[a][g][p][e],
              dist_ll_num_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          dist_rr_num[a][g][p][e] = new TH1D(dist_rr_num_n[a][g][p][e],
              dist_rr_num_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          dist_ll_den[a][g][p][e] = new TH1D(dist_ll_den_n[a][g][p][e],
              dist_ll_den_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          dist_rr_den[a][g][p][e] = new TH1D(dist_rr_den_n[a][g][p][e],
              dist_rr_den_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          
          dist_ll_num_e[a][g][p][e] = new TH1D(dist_ll_num_e_n[a][g][p][e],
              dist_ll_num_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          dist_rr_num_e[a][g][p][e] = new TH1D(dist_rr_num_e_n[a][g][p][e],
              dist_rr_num_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          dist_ll_den_e[a][g][p][e] = new TH1D(dist_ll_den_e_n[a][g][p][e],
              dist_ll_den_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          dist_rr_den_e[a][g][p][e] = new TH1D(dist_rr_den_e_n[a][g][p][e],
              dist_rr_den_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);

          numer[a][g][p][e] = new TH1D(numer_n[a][g][p][e],
              numer_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          denom[a][g][p][e] = new TH1D(denom_n[a][g][p][e],
              denom_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          numer_e[a][g][p][e] = new TH1D(numer_e_n[a][g][p][e],
              numer_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          denom_e[a][g][p][e] = new TH1D(denom_e_n[a][g][p][e],
              denom_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);

          numer_e_sqrt[a][g][p][e] = new TH1D(numer_e_sqrt_n[a][g][p][e],
              numer_e_sqrt_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);

          asym[a][g][p][e] = new TH1D(asym_n[a][g][p][e],
              asym_t[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);
          asym_e[a][g][p][e] = new TH1D(asym_e_n[a][g][p][e],
              asym_e_n[a][g][p][e],phi_bins,RD->env->PhiLow,RD->env->PhiHigh);

          // build left & right terms
          if(a==3)
          {
            dist_ll_num[a][g][p][e]->Add(phi_dist_num[a][0][g][p][e],phi_dist_num[a][3][g][p][e],1.0,1.0);
            dist_rr_num[a][g][p][e]->Add(phi_dist_num[a][1][g][p][e],phi_dist_num[a][2][g][p][e],1.0,1.0);
            dist_ll_den[a][g][p][e]->Add(phi_dist_den[a][0][g][p][e],phi_dist_den[a][3][g][p][e],1.0,1.0);
            dist_rr_den[a][g][p][e]->Add(phi_dist_den[a][1][g][p][e],phi_dist_den[a][2][g][p][e],1.0,1.0);
            dist_ll_num_e[a][g][p][e]->Add(phi_dist_num_e[a][0][g][p][e],phi_dist_num_e[a][3][g][p][e],1.0,1.0);
            dist_rr_num_e[a][g][p][e]->Add(phi_dist_num_e[a][1][g][p][e],phi_dist_num_e[a][2][g][p][e],1.0,1.0);
            dist_ll_den_e[a][g][p][e]->Add(phi_dist_den_e[a][0][g][p][e],phi_dist_den_e[a][3][g][p][e],1.0,1.0);
            dist_rr_den_e[a][g][p][e]->Add(phi_dist_den_e[a][1][g][p][e],phi_dist_den_e[a][2][g][p][e],1.0,1.0);
          }
          else if(a==1)
          {
            dist_ll_num[a][g][p][e]->Add(phi_dist_num[a][1][g][p][e],phi_dist_num[a][3][g][p][e],1.0,1.0);
            dist_rr_num[a][g][p][e]->Add(phi_dist_num[a][0][g][p][e],phi_dist_num[a][2][g][p][e],1.0,1.0);
            dist_ll_den[a][g][p][e]->Add(phi_dist_den[a][1][g][p][e],phi_dist_den[a][3][g][p][e],1.0,1.0);
            dist_rr_den[a][g][p][e]->Add(phi_dist_den[a][0][g][p][e],phi_dist_den[a][2][g][p][e],1.0,1.0);
            dist_ll_num_e[a][g][p][e]->Add(phi_dist_num_e[a][1][g][p][e],phi_dist_num_e[a][3][g][p][e],1.0,1.0);
            dist_rr_num_e[a][g][p][e]->Add(phi_dist_num_e[a][0][g][p][e],phi_dist_num_e[a][2][g][p][e],1.0,1.0);
            dist_ll_den_e[a][g][p][e]->Add(phi_dist_den_e[a][1][g][p][e],phi_dist_den_e[a][3][g][p][e],1.0,1.0);
            dist_rr_den_e[a][g][p][e]->Add(phi_dist_den_e[a][0][g][p][e],phi_dist_den_e[a][2][g][p][e],1.0,1.0);
          }
          else if(a==2)
          {
            dist_ll_num[a][g][p][e]->Add(phi_dist_num[a][2][g][p][e],phi_dist_num[a][3][g][p][e],1.0,1.0);
            dist_rr_num[a][g][p][e]->Add(phi_dist_num[a][0][g][p][e],phi_dist_num[a][1][g][p][e],1.0,1.0);
            dist_ll_den[a][g][p][e]->Add(phi_dist_den[a][2][g][p][e],phi_dist_den[a][3][g][p][e],1.0,1.0);
            dist_rr_den[a][g][p][e]->Add(phi_dist_den[a][0][g][p][e],phi_dist_den[a][1][g][p][e],1.0,1.0);
            dist_ll_num_e[a][g][p][e]->Add(phi_dist_num_e[a][2][g][p][e],phi_dist_num_e[a][3][g][p][e],1.0,1.0);
            dist_rr_num_e[a][g][p][e]->Add(phi_dist_num_e[a][0][g][p][e],phi_dist_num_e[a][1][g][p][e],1.0,1.0);
            dist_ll_den_e[a][g][p][e]->Add(phi_dist_den_e[a][2][g][p][e],phi_dist_den_e[a][3][g][p][e],1.0,1.0);
            dist_rr_den_e[a][g][p][e]->Add(phi_dist_den_e[a][0][g][p][e],phi_dist_den_e[a][1][g][p][e],1.0,1.0);
          }

          // build numer and denom
          numer[a][g][p][e]->Add(dist_ll_num[a][g][p][e],dist_rr_num[a][g][p][e],1.0,-1.0);
          denom[a][g][p][e]->Add(dist_ll_den[a][g][p][e],dist_rr_den[a][g][p][e],1.0,1.0);
          numer_e[a][g][p][e]->Add(dist_ll_num_e[a][g][p][e],dist_rr_num_e[a][g][p][e],1.0,1.0); // n.b.: ll & rr summed
          denom_e[a][g][p][e]->Add(dist_ll_den_e[a][g][p][e],dist_rr_den_e[a][g][p][e],1.0,1.0);
          
          // numer_e_sqrt = sqrt(numer_e)
          for(Int_t b=1; b<=numer_e[a][g][p][e]->GetNbinsX(); b++)
          {
            bc_e = numer_e[a][g][p][e]->GetBinContent(b);
            bc_e = sqrt(bc_e);
            numer_e_sqrt[a][g][p][e]->SetBinContent(b,bc_e);
          };

          asym[a][g][p][e]->Divide(numer[a][g][p][e],denom[a][g][p][e],1.0,1.0);
          asym_e[a][g][p][e]->Divide(numer_e_sqrt[a][g][p][e],denom_e[a][g][p][e],1.0,1.0);

          printf(asym[a][g][p][e]->GetTitle());
          printf("\n");

          
          // set statistical errors of asym bins using values computed in asym_e
          for(Int_t b=1; b<=asym[a][g][p][e]->GetNbinsX(); b++)
          {
            bc_e = asym_e[a][g][p][e]->GetBinContent(b);
            // compare computed stat error with ROOT's error propagation
            //bc_ee = asym[a][g][p][e]->GetBinError(b);
            //gSystem->RedirectOutput("error_check.dat","a");
            //printf("%.20f\n",bc_ee-bc_e); 
            //gSystem->RedirectOutput(0);
            asym[a][g][p][e]->SetBinError(b,bc_e);
          };

          // n.b. for one phi bin, constant fit & error matches the bin & its error
          sprintf(asym_fit_n[a][g][p][e],"asym_fit_a%d_g%d_p%d_e%d",a,g,p,e);
          if(a==3) 
          {
            if(ANALYSIS_TYPE==kTrans)
            {
              asym_fit[a][g][p][e] = new TF1(asym_fit_n[a][g][p][e],"[0]+[1]*cos(2*x)",RD->env->PhiLow,RD->env->PhiHigh);
              asym_fit[a][g][p][e]->SetParName(0,"A_{#Sigma}");
              asym_fit[a][g][p][e]->SetParName(1,"A_{TT}");
            }
            else if(ANALYSIS_TYPE==kLong)
            {
              asym_fit[a][g][p][e] = new TF1(asym_fit_n[a][g][p][e],"pol0",RD->env->PhiLow,RD->env->PhiHigh);
              asym_fit[a][g][p][e]->SetParName(0,"A_{LL}");
            };
          }
          else if(a==1 || a==2)
          {
            if(ANALYSIS_TYPE==kTrans)
            {
              asym_fit[a][g][p][e] = new TF1(asym_fit_n[a][g][p][e],"[0]+[1]*cos(x)",RD->env->PhiLow,RD->env->PhiHigh);
              asym_fit[a][g][p][e]->SetParName(0,"R");
              if(a==1) asym_fit[a][g][p][e]->SetParName(1,"A_{N}^{Y}");
              else if(a==2) asym_fit[a][g][p][e]->SetParName(1,"A_{N}^{B}");
            }
            else if(ANALYSIS_TYPE==kLong)
            {
              asym_fit[a][g][p][e] = new TF1(asym_fit_n[a][g][p][e],"pol0",RD->env->PhiLow,RD->env->PhiHigh);
              if(a==1) asym_fit[a][g][p][e]->SetParName(0,"A_{L}^{Y}");
              else if(a==2) asym_fit[a][g][p][e]->SetParName(0,"A_{L}^{B}");
            };
          };
          asym[a][g][p][e]->Fit(asym_fit[a][g][p][e],"Q","",RD->env->PhiLow,RD->env->PhiHigh);
          if(asym_fit[a][g][p][e]!=NULL)
          {
            if(ANALYSIS_TYPE==kTrans) qq=1;
            else if(ANALYSIS_TYPE==kLong) qq=0;
            if(a==3) asym_value[a][g][p][e]+=asym_fit[a][g][p][e]->GetParameter(qq);
            else if(a==1 || a==2) asym_value[a][g][p][e]+=asym_fit[a][g][p][e]->GetParameter(qq);
            asym_value_cnt[a][g][p][e]++;
          };
        };
      };
    };
  };


  // set plot ranges
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          //asym[a][g][p][e]->GetYaxis()->SetRangeUser(2*asym_min[a][g][p][e],2*asym_max[a][g][p][e]);
          asym[a][g][p][e]->GetXaxis()->SetRangeUser(RD->env->PhiLow,RD->env->PhiHigh);
        };
      };
    };
  };


  // average asym_value for each bin
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          asym_value[a][g][p][e] /= ((Float_t)asym_value_cnt[a][g][p][e]);
          printf("g%d p%d e%d <%s>=%f\n",g,p,e,asym_title[a],asym_value[a][g][p][e]);
        };
      };
    };
  };


  // open wdists
  TH1D * pt_wdist[eta_bins][en_bins]; // [eta] [en]
  TH1D * en_wdist[eta_bins][pt_bins]; // [eta] [en]
  char pt_wdist_n[eta_bins][en_bins][64];
  char en_wdist_n[eta_bins][pt_bins][64];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      sprintf(pt_wdist_n[g][e],"pt_wdist_tot_%s_g%d_e%d",evclass,g,e);
      pt_wdist[g][e] = (TH1D*) infile->Get(pt_wdist_n[g][e]);
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      sprintf(en_wdist_n[g][p],"en_wdist_tot_%s_g%d_p%d",evclass,g,p);
      en_wdist[g][p] = (TH1D*) infile->Get(en_wdist_n[g][p]);
    };
  };


  // wdist weighting histograms
  Float_t pt_cc[pt_bins];
  Float_t pt_ww[pt_bins];
  Float_t en_cc[pt_bins];
  Float_t en_ww[pt_bins];
  Int_t NWBINS = pt_wdist[0][0]->GetNbinsX();
  TH1D * pt_wdist_sub[eta_bins][en_bins][pt_bins];
  TH1D * en_wdist_sub[eta_bins][pt_bins][en_bins];
  char pt_wdist_sub_n[eta_bins][en_bins][pt_bins][64];
  char en_wdist_sub_n[eta_bins][pt_bins][en_bins][64];
  Int_t iter = 0;
  // first create the wdist "sub" distributions
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        sprintf(pt_wdist_sub_n[g][e][p],"pt_wdist_sub_g%d_e%d_ptbin%d",g,e,p);
        pt_wdist_sub[g][e][p] = new TH1D(pt_wdist_sub_n[g][e][p],
          pt_wdist_sub_n[g][e][p],NWBINS,RD->env->PtLow,RD->env->PtHigh);
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        sprintf(en_wdist_sub_n[g][p][e],"en_wdist_sub_g%d_p%d_enbin%d",g,p,e);
        en_wdist_sub[g][p][e] = new TH1D(en_wdist_sub_n[g][p][e],
          en_wdist_sub_n[g][p][e],NWBINS,RD->env->EnLow,RD->env->EnHigh);
      };
    };
  };
  // then fill the wdist "sub" distributions
  Double_t bincont,bincent;
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t b=1; b<=NWBINS; b++)
      {
        bincont = pt_wdist[g][e]->GetBinContent(b);
        bincent = pt_wdist[g][e]->GetBinCenter(b);
        for(Int_t p=0; p<pt_bins; p++)
        {
          if(bincent>=RD->env->PtDiv(p) && bincent<RD->env->PtDiv(p+1))
          {
            pt_wdist_sub[g][e][p]->SetBinContent(b,bincont);
          };
        };
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t b=1; b<=NWBINS; b++)
      {
        bincont = en_wdist[g][p]->GetBinContent(b);
        bincent = en_wdist[g][p]->GetBinCenter(b);
        for(Int_t e=0; e<en_bins; e++)
        {
          if(bincent>=RD->env->EnDiv(e) && bincent<RD->env->EnDiv(e+1))
          {
            en_wdist_sub[g][p][e]->SetBinContent(b,bincont);
          };
        };
      };
    };
  };



  // draw TLines to wdists which indicate the weighting
  TLine * pt_div_line[eta_bins][en_bins][pt_bins+1];
  TLine * en_div_line[eta_bins][pt_bins][en_bins+1];
  TLine * pt_cent_line[eta_bins][en_bins][pt_bins];
  TLine * en_cent_line[eta_bins][pt_bins][en_bins];
  TLine * pt_width_line[eta_bins][en_bins][pt_bins];
  TLine * en_width_line[eta_bins][pt_bins][en_bins];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        pt_div_line[g][e][p] = new TLine(RD->env->PtDiv(p),0,RD->env->PtDiv(p),pt_wdist[g][e]->GetMaximum());
        pt_cent_line[g][e][p] = new TLine(pt_wdist_sub[g][e][p]->GetMean(),
                                          pt_wdist[g][e]->GetMaximum() * 1/4,
                                          pt_wdist_sub[g][e][p]->GetMean(),
                                          pt_wdist[g][e]->GetMaximum() * 3/4);
        pt_width_line[g][e][p] = new TLine(pt_wdist_sub[g][e][p]->GetMean() - pt_wdist_sub[g][e][p]->GetRMS(),
                                           pt_wdist[g][e]->GetMaximum()/2,
                                           pt_wdist_sub[g][e][p]->GetMean() + pt_wdist_sub[g][e][p]->GetRMS(),
                                           pt_wdist[g][e]->GetMaximum()/2);
        pt_div_line[g][e][p]->SetLineColor(kBlack);
        pt_cent_line[g][e][p]->SetLineColor(kRed);
        pt_width_line[g][e][p]->SetLineColor(kRed);
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        en_div_line[g][p][e] = new TLine(RD->env->EnDiv(e),0,RD->env->EnDiv(e),en_wdist[g][p]->GetMaximum());
        en_cent_line[g][p][e] = new TLine(en_wdist_sub[g][p][e]->GetMean(),
                                          en_wdist[g][p]->GetMaximum() * 1/4,
                                          en_wdist_sub[g][p][e]->GetMean(),
                                          en_wdist[g][p]->GetMaximum() * 3/4);
        en_width_line[g][p][e] = new TLine(en_wdist_sub[g][p][e]->GetMean() - en_wdist_sub[g][p][e]->GetRMS(),
                                           en_wdist[g][p]->GetMaximum()/2,
                                           en_wdist_sub[g][p][e]->GetMean() + en_wdist_sub[g][p][e]->GetRMS(),
                                           en_wdist[g][p]->GetMaximum()/2);
        en_div_line[g][p][e]->SetLineColor(kBlack);
        en_cent_line[g][p][e]->SetLineColor(kRed);
        en_width_line[g][p][e]->SetLineColor(kRed);
      };
    };
  };
  

  // make TCanvases for wdist plots with marker lines
  TCanvas * pt_wdist_canv[eta_bins][en_bins];
  TCanvas * en_wdist_canv[eta_bins][pt_bins];
  char pt_wdist_canv_n[eta_bins][en_bins][64];
  char en_wdist_canv_n[eta_bins][pt_bins][64];
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      sprintf(pt_wdist_canv_n[g][e],"pt_wdist_canv_g%d_e%d",g,e);
      pt_wdist_canv[g][e] = new TCanvas(pt_wdist_canv_n[g][e],pt_wdist_canv_n[g][e],700,500);
      pt_wdist[g][e]->Draw();
      for(Int_t p=0; p<pt_bins; p++)
      {
        pt_div_line[g][e][p]->Draw();
        pt_cent_line[g][e][p]->Draw();
        pt_width_line[g][e][p]->Draw();
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      sprintf(en_wdist_canv_n[g][p],"en_wdist_canv_g%d_p%d",g,p);
      en_wdist_canv[g][p] = new TCanvas(en_wdist_canv_n[g][p],en_wdist_canv_n[g][p],700,500);
      en_wdist[g][p]->Draw();
      for(Int_t e=0; e<en_bins; e++)
      {
        en_div_line[g][p][e]->Draw();
        en_cent_line[g][p][e]->Draw();
        en_width_line[g][p][e]->Draw();
      };
    };
  };



  // kinematic dependence plots
  TGraphErrors * en_dep[NPARAM][asym_bins][eta_bins][pt_bins]; // en dependent plots, one for each pt bin (statistical errors)
  TGraphErrors * pt_dep[NPARAM][asym_bins][eta_bins][en_bins]; // pt dependent plots, one for each en bin
  char en_dep_t[NPARAM][asym_bins][eta_bins][pt_bins][256];
  char pt_dep_t[NPARAM][asym_bins][eta_bins][en_bins][256];
  Int_t en_dep_cnt[NPARAM][asym_bins][eta_bins][pt_bins]; // point counter
  Int_t pt_dep_cnt[NPARAM][asym_bins][eta_bins][en_bins];
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          en_dep_cnt[z][a][g][p]=0;
        };
        for(Int_t e=0; e<en_bins; e++) 
        {
          pt_dep_cnt[z][a][g][e]=0;
        };
      };
    };
  };

  Double_t val_en[NPARAM][asym_bins][eta_bins][pt_bins][en_bins];     // arrays for en dependent plots, one for each pt bin
  Double_t err_en[NPARAM][asym_bins][eta_bins][pt_bins][en_bins];     // (statistical error)
  Double_t cent_en[NPARAM][asym_bins][eta_bins][pt_bins][en_bins];    // (energy bin center)
  Double_t width_en[NPARAM][asym_bins][eta_bins][pt_bins][en_bins];   // (energy bin width)
  Double_t zeroz_en[NPARAM][asym_bins][eta_bins][pt_bins][en_bins];   // array of zeros

  Double_t val_pt[NPARAM][asym_bins][eta_bins][en_bins][pt_bins];     // arrays for pt dependent plots, one
  Double_t err_pt[NPARAM][asym_bins][eta_bins][en_bins][pt_bins];     // (statistical error)
  Double_t cent_pt[NPARAM][asym_bins][eta_bins][en_bins][pt_bins];    // (pt bin center)
  Double_t width_pt[NPARAM][asym_bins][eta_bins][en_bins][pt_bins];   // (pt bin width)
  Double_t zeroz_pt[NPARAM][asym_bins][eta_bins][en_bins][pt_bins];   // array of zeros

  // en dependent points for each pt bin
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            zeroz_en[z][a][g][p][e]=0;
            if(asym_fit[a][g][p][e])
            {
              val_en[z][a][g][p][en_dep_cnt[z][a][g][p]] = asym_fit[a][g][p][e]->GetParameter(z);
              // error of fit to asym; asym statistical errors set via asym_e
                err_en[z][a][g][p][en_dep_cnt[z][a][g][p]] = asym_fit[a][g][p][e]->GetParError(z);
              // estimated statistical error
                //err_en[z][a][g][p][en_dep_cnt[z][a][g][p]]=1/(0.55*0.55)*1/sqrt(yield[0][g][p][e]+yield[1][g][p][e]+yield[2][g][p][e]+yield[3][g][p][e]);

              //cent_en[z][a][g][p][en_dep_cnt[z][a][g][p]] = RD->env->EnDiv(e) + ((RD->env->EnDiv(e+1)-RD->env->EnDiv(e))/2.0);
              //width_en[z][a][g][p][en_dep_cnt[z][a][g][p]] = (RD->env->EnDiv(e+1)-RD->env->EnDiv(e))/2.0;
              cent_en[z][a][g][p][en_dep_cnt[z][a][g][p]] = en_wdist_sub[g][p][e]->GetMean();
              width_en[z][a][g][p][en_dep_cnt[z][a][g][p]] = en_wdist_sub[g][p][e]->GetRMS();

              en_dep_cnt[z][a][g][p]++;
            };
          };
          // asymmetry vs. en with statistical error bars
            en_dep[z][a][g][p] = new TGraphErrors(en_dep_cnt[z][a][g][p],cent_en[z][a][g][p],
              val_en[z][a][g][p],width_en[z][a][g][p],err_en[z][a][g][p]);
            sprintf(en_dep_t[z][a][g][p],
              "%s #pm #sigma %s vs. E for p_{T}#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",
              asym_title_kd[z][a],asym_title_kd[z][a],
              RD->env->PtDiv(p),RD->env->PtDiv(p+1),RD->env->EtaDiv(g),RD->env->EtaDiv(g+1));
          
          en_dep[z][a][g][p]->SetTitle(en_dep_t[z][a][g][p]);
          en_dep[z][a][g][p]->GetXaxis()->SetTitle("E (GeV)");
        };
      };
    };
  };

  // pt dependent points for each en bin
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          for(Int_t p=0; p<pt_bins; p++)
          {
            zeroz_pt[z][a][g][e][p]=0;
            if(asym_fit[a][g][p][e])
            {
              val_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]] = asym_fit[a][g][p][e]->GetParameter(z);
              // error of fit to asym; asym statistical errors set via asym_e
                err_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]] = asym_fit[a][g][p][e]->GetParError(z);
              // estimated statistical error
                //err_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]]=1/(0.55*0.55)*1/sqrt(yield[0][g][p][e]+yield[1][g][p][e]+yield[2][g][p][e]+yield[3][g][p][e]);

              //cent_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]] = RD->env->PtDiv(p) + ((RD->env->PtDiv(p+1)-RD->env->PtDiv(p))/2.0);
              //width_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]] = (RD->env->PtDiv(p+1)-RD->env->PtDiv(p))/2.0;
              cent_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]] = pt_wdist_sub[g][e][p]->GetMean();
              width_pt[z][a][g][e][pt_dep_cnt[z][a][g][e]] = pt_wdist_sub[g][e][p]->GetRMS();

              pt_dep_cnt[z][a][g][e]++;
            };
          };
          // asymmetry vs. pt with statistical error bars
            pt_dep[z][a][g][e] = new TGraphErrors(pt_dep_cnt[z][a][g][e],cent_pt[z][a][g][e],
              val_pt[z][a][g][e],width_pt[z][a][g][e],err_pt[z][a][g][e]);
            sprintf(pt_dep_t[z][a][g][e],"%s #pm #sigma %s vs. p_{T} for E#in[%.2f,%.2f) and #eta#in[%.2f,%.2f)",
              asym_title_kd[z][a],asym_title_kd[z][a],
              RD->env->EnDiv(e),RD->env->EnDiv(e+1),RD->env->EtaDiv(g),RD->env->EtaDiv(g+1));

          pt_dep[z][a][g][e]->SetTitle(pt_dep_t[z][a][g][e]);
          pt_dep[z][a][g][e]->GetXaxis()->SetTitle("p_{T} (GeV)");
        };
      };
    };
  };



  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          en_dep[z][a][g][p]->GetYaxis()->SetTitle(asym_title_kd[z][a]);
          en_dep[z][a][g][p]->SetMarkerStyle(kFullCircle);
          en_dep[z][a][g][p]->SetMarkerColor(kRed);
          en_dep[z][a][g][p]->GetYaxis()->SetTitleOffset(1.5);
        };
        for(Int_t e=0; e<en_bins; e++)
        {
          pt_dep[z][a][g][e]->GetYaxis()->SetTitle(asym_title_kd[z][a]);
          pt_dep[z][a][g][e]->SetMarkerStyle(kFullCircle);
          pt_dep[z][a][g][e]->SetMarkerColor(kRed);
          pt_dep[z][a][g][e]->GetYaxis()->SetTitleOffset(1.5);
        };
      };
    };
  };


  // write phi dists
  printf("writing spin.root...\n");
  TFile * outfile = new TFile("spin.root","RECREATE");
  if(ANALYSIS_TYPE==kTrans) {
    outfile->mkdir("A_Sigma");
    outfile->mkdir("R_blue");
    outfile->mkdir("R_yellow");
    outfile->mkdir("A_TT");
    outfile->mkdir("A_N_blue");
    outfile->mkdir("A_N_yellow");
  }
  else if(ANALYSIS_TYPE==kLong) {
    outfile->mkdir("A_LL");
    outfile->mkdir("A_L_blue");
    outfile->mkdir("A_L_yellow");
  };
  char en_dep_n[NPARAM][asym_bins][eta_bins][pt_bins][32];
  char pt_dep_n[NPARAM][asym_bins][eta_bins][en_bins][32];
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++)
    {
      if(z==0)
      {
        if(ANALYSIS_TYPE==kTrans) {
          if(a==3) outfile->cd("/A_Sigma");
          else if(a==1) outfile->cd("/R_yellow");
          else if(a==2) outfile->cd("/R_blue");
        }
        else if(ANALYSIS_TYPE==kLong) {
          if(a==3) outfile->cd("/A_LL");
          else if(a==1) outfile->cd("/A_L_yellow");
          else if(a==2) outfile->cd("/A_L_blue");
        }
      }
      else if(z==1)
      {
        if(a==3) outfile->cd("/A_TT");
        else if(a==1) outfile->cd("/A_N_yellow");
        else if(a==2) outfile->cd("/A_N_blue");
      }:
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          sprintf(en_dep_n[z][a][g][p],"en_dep_z%d_a%d_g%d_p%d",z,a,g,p);
          en_dep[z][a][g][p]->Write(en_dep_n[z][a][g][p]);
        };
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(pt_dep_n[z][a][g][e],"pt_dep_z%d_a%d_g%d_e%d",z,a,g,e);
          pt_dep[z][a][g][e]->Write(pt_dep_n[z][a][g][e]);
        };
      };
    };
  };

  // write phi_dist_num/den and asym to both shift and amplitude asymmetry directories
  for(Int_t a=1; a<asym_bins; a++)
  {
    for(z=0; z<NPARAM; z++) {
      if(z==0) {
        if(ANALYSIS_TYPE==kTrans) {
          if(a==3) outfile->cd("/A_Sigma");
          else if(a==1) outfile->cd("/R_yellow");
          else if(a==2) outfile->cd("/R_blue");
        }
        else if(ANALYSIS_TYPE==kLong) {
          if(a==3) outfile->cd("/A_LL");
          else if(a==1) outfile->cd("/A_L_yellow");
          else if(a==2) outfile->cd("/A_L_blue");
        }
      }
      else if(z==1) {
        if(a==3) outfile->cd("/A_TT");
        else if(a==1) outfile->cd("/A_N_yellow");
        else if(a==2) outfile->cd("/A_N_blue");
      };
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            for(Int_t s=0; s<4; s++)
            {
              phi_dist_num[a][s][g][p][e]->Write(phi_dist_num_n[a][s][g][p][e]);
              phi_dist_den[a][s][g][p][e]->Write(phi_dist_den_n[a][s][g][p][e]);
            };
          };
        };
      };
    };
  };

  for(Int_t a=1; a<asym_bins; a++)
  {
    for(z=0; z<NPARAM; z++) {
      if(z==0) {
        if(ANALYSIS_TYPE==kTrans) {
          if(a==3) outfile->cd("/A_Sigma");
          else if(a==1) outfile->cd("/R_yellow");
          else if(a==2) outfile->cd("/R_blue");
        }
        else if(ANALYSIS_TYPE==kLong) {
          if(a==3) outfile->cd("/A_LL");
          else if(a==1) outfile->cd("/A_L_yellow");
          else if(a==2) outfile->cd("/A_L_blue");
        }
      }
      else if(z==1) {
        if(a==3) outfile->cd("/A_TT");
        else if(a==1) outfile->cd("/A_N_yellow");
        else if(a==2) outfile->cd("/A_N_blue");
      };

      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            asym[a][g][p][e]->Write(asym_n[a][g][p][e]);
          };
        };
      };
    };
  };
  /*
  for(Int_t a=1; a<asym_bins; a++)
  {
    if(a==3) outfile->cd("/A_Sigma");
    else if(a==1) outfile->cd("/R_yellow");
    else if(a==2) outfile->cd("/R_blue");
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          for(Int_t s=0; s<4; s++)
          {
            phi_dist_num[a][s][g][p][e]->Write(phi_dist_num_n[a][s][g][p][e]);
            phi_dist_den[a][s][g][p][e]->Write(phi_dist_den_n[a][s][g][p][e]);
          };
        };
      };
    };
  };
  for(Int_t a=1; a<asym_bins; a++)
  {
    if(a==3) outfile->cd("/A_Sigma");
    else if(a==1) outfile->cd("/R_yellow");
    else if(a==2) outfile->cd("/R_blue");
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          asym[a][g][p][e]->Write(asym_n[a][g][p][e]);
        };
      };
    };
  };
  */
  

  // write out wdist information in directory "bin_weighting"
  outfile->mkdir("bin_weighting");
  outfile->cd("/bin_weighting");
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++) pt_wdist_canv[g][e]->Write();
    for(Int_t p=0; p<pt_bins; p++) en_wdist_canv[g][p]->Write();
  };
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++) pt_wdist[g][e]->Write();
    for(Int_t p=0; p<pt_bins; p++) en_wdist[g][p]->Write();
  }
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t p=0; p<pt_bins; p++) pt_wdist_sub[g][e][p]->Write();
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++) en_wdist_sub[g][p][e]->Write();
    };
  };

  printf("written\n");
};
