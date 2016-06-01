#ifndef BBCtiles_
#define BBCtiles_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include "TSystem.h"
#include "TMath.h"
#include "TH2Poly.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMatrix.h"
#include "TVector.h"


class BBCtiles : public TObject
{
  public:
    BBCtiles();
    TH2Poly * MakeNewPoly(Int_t sl0);
      
    void HexToCart(Double_t a0, Int_t x0, Int_t y0, Int_t z0, Double_t &xx0, Double_t &yy0);

    Int_t GetBinOfTile(Int_t tile0);
    Int_t GetSlOfTile(Int_t tile0);
    Int_t GetSlOfPMT(Int_t pmt0);
    Int_t GetPMTOfTile(Int_t tile0);
    Int_t GetTileOfPMT(Int_t pmt0, Int_t whichTile);
    Int_t GetXhexOfTile(Int_t tile0);
    Int_t GetYhexOfTile(Int_t tile0);
    Int_t GetZhexOfTile(Int_t tile0);

    Double_t GetAzimuthOfTile(Int_t tile0,Bool_t returnDeg=0);
    Double_t GetAveAzimuthOfPMT(Int_t pmt0);

    Double_t ComputeAzimuthOfTile(Int_t tile0,Bool_t returnDeg=0);
    Double_t ComputeAveAzimuthOfPMT(Int_t pmt0);

    void InitRecenteringValues();

    void DrawBBC();
    void PrintBBC();

    void UpdateEvent();
    void ResetEvent();

    void ComputeEVP(Int_t ew0, Int_t sl0); // sets EVP, Xfow, and Yflow
    void ComputeMoments(Int_t ew0, Int_t sl0);

    void DrawEvent();

    Bool_t IsVertical();
    Bool_t IsHorizontal();


    TH2Poly * tile_poly[2]; // [sl (sl=0 for small, sl=1 for large)]
    TH2Poly * pmt_poly[2]; // [sl]
    TH2Poly * adc_poly[2][2]; // [ew] [sl] // single event ADC
    TH2Poly * tac_poly[2][2]; // [ew] [sl] // single event TAC
    TH2Poly * acc_ev_poly[2][2]; // [ew] [sl] // accumulated hit display


    /* event variables */
    Char_t QTN[2][2]; // [ew] [sl]
    Char_t Idx[2][2][16]; // [ew] [sl] [channel]
    Short_t ADC[2][2][16]; // [ew] [sl] [channel]
    Short_t TAC[2][2][16]; // [ew] [sl] [channel]
    Float_t vertex;

    
    Double_t dir[2][2]; // [ew] [sl] // first fourier coeff
    Double_t EVP[2][2]; // [ew] [sl] // second fourier coeff = event plane azimuth
    Double_t Xflow[2][2]; // [ew] [sl]
    Double_t Yflow[2][2]; // [ew] [sl]

    Double_t Xflow_cor[2][2]; // [ew] [sl] // re-centered
    Double_t Yflow_cor[2][2]; // [ew] [sl] // re-centered
    Double_t EVP_cor[2][2]; // [ew] [sl] // re-centered

    // momoent variables for planarity
    Double_t xbar[2][2]; // [ew] [sl];
    Double_t ybar[2][2];
    Double_t esum[2][2];
    Double_t sigma_x[2][2];
    Double_t sigma_y[2][2];
    Double_t sigma_xy[2][2];
    Double_t sigma_min[2][2];
    Double_t sigma_max[2][2];
    Double_t discrim;
    Double_t sigma_theta[2][2]; // angle of sigma_max eigenvector


    // canvases
    TCanvas * bbc_canv;
    TCanvas * ev_canv;


    
  protected:
    int ew,sl,tt,pp,tpxyz;
    Bool_t initTlines;
    Double_t xfl1,yfl1,xfl2,yfl2;

    // BBC coordinates map
    // [0=small 1=large] [tile enumerator] [0=tile# 1=pmt# 2=xhex 3=yhex 4=zhex]
    // n.b. TH2Poly bin# = tile# - 18 * sl (where sl=0 for small, 1 for large)
    // - see below for hex coordinates
    Int_t coord[2][18][5];


    /* hexagonal coordinates
     *  
     *      Y
     *     
     *        o 
     *         \     /
     *          \   /
     *           \ /
     *        ----o----o  X
     *           / \
     *          /   \
     *         /     \
     *        o
     *     
     *      Z
     *
     */


    // fast mapping of pmt number to tiles; since there are often more than
    // 1 tile for a given pmt number (up to 3 tiles), we save all tile numbers
    // for each pmt (or zero if nonexistant)
    Int_t coord_pmt[25][3]; // [pmt number] [tile1, tile2 tile3];


    // tile or average pmt azimuths (for fast accessors)
    Double_t azi_tile[37]; // [tile number]
    Double_t azi_pmt[25]; // [pmt number]


    Double_t cellsize[2]; // [sl] // event display hexagon side size
    TLine * hexline[6][2][18]; // [hex side] [sl] [tile]
    Int_t hexlinecolor[2]; // [sl]

    Double_t xcart,ycart; // cartesian coordinates for event display
    Int_t xhex,yhex,zhex; // hex coordinates
    Double_t x[6];
    Double_t y[6];

    TMatrix * sigma_mtx;
    TVector * sigma_evals;
    Double_t sigma_min_tmp,sigma_max_tmp;

    // EVP directionality
    Double_t left_sum,right_sum;
    Int_t directionality[2][2]; // 0 for positive, 1 for negative, -1 for undetermined, 2 for left/right equal
    Double_t left_sum_cor,right_sum_cor; // (for re-centered evp calculation)
    Int_t directionality_cor[2][2]; // (for re-centered evp calculation)


    Bool_t PMTmasked[25]; // [pmt] // if true, don't use for EVP calculation

    TLine * evp_line[2][2]; // [ew] [sl]
    TLine * evp_cor_line[2][2]; // [ew] [sl]
    Double_t scale;

    TString datastring[2][2]; // [ew] [sl]

    Double_t Xflow_calib[2][2][16][2]; // [ew] [sl] [qtn-1] [0=mean, 1=rms]
    Double_t Yflow_calib[2][2][16][2]; // [ew] [sl] [qtn-1] [0=mean, 1=rms]

    Double_t div[4];
    Double_t eee[2];



    ClassDef(BBCtiles,1);

};
#endif
