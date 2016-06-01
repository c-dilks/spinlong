#include "BBCtiles.h"

namespace
{
  enum ew_enum {kE,kW};
  enum sl_enum {kS,kL};
  enum tpxyz_enum {kTile,kPMT,kX,kY,kZ};
  const Double_t sq = TMath::Sqrt(3);
};


BBCtiles::BBCtiles()
{
  // BBC coordinate map
  // (first define it here as constructor-local array)
  Int_t coord0[2][18][5] = 
  {
    {
      {  1,  1,  0,  1, -1 },
      {  2,  2,  1,  0, -1 },
      {  3,  3,  1, -1,  0 },
      {  4,  4,  0, -1,  1 },
      {  5,  5, -1,  0,  1 },
      {  6,  6, -1,  1,  0 },
      {  7,  7, -1,  2, -1 },
      {  8,  8,  0,  2, -2 },
      {  9,  7,  1,  1, -2 },
      { 10,  9,  2,  0, -2 },
      { 11, 10,  2, -1, -1 },
      { 12, 11,  2, -2,  0 },
      { 13, 12,  1, -2,  1 },
      { 14, 13,  0, -2,  2 },
      { 15, 12, -1, -1,  2 },
      { 16, 14, -2,  0,  2 },
      { 17, 15, -2,  1,  1 },
      { 18, 16, -2,  2,  0 }
    },
    {
      { 19, 17,  0,  1, -1 },
      { 20, 18,  1,  0, -1 },
      { 21, 18,  1, -1,  0 },
      { 22, 19,  0, -1,  1 },
      { 23, 20, -1,  0,  1 },
      { 24, 20, -1,  1,  0 },
      { 25, 21, -1,  2, -1 },
      { 26, 21,  0,  2, -2 },
      { 27, 21,  1,  1, -2 },
      { 28, 22,  2,  0, -2 },
      { 29, 22,  2, -1, -1 },
      { 30, 22,  2, -2,  0 },
      { 31, 23,  1, -2,  1 },
      { 32, 23,  0, -2,  2 },
      { 33, 23, -1, -1,  2 },
      { 34, 24, -2,  0,  2 },
      { 35, 24, -2,  1,  1 },
      { 36, 24, -2,  2,  0 }
    }
  };


  // copy constructor-local array to class member "coord"
  for(sl=0; sl<2; sl++)
  {
    for(tt=0; tt<18; tt++)
    {
      for(tpxyz=0; tpxyz<5; tpxyz++)
      {
        coord[sl][tt][tpxyz] = coord0[sl][tt][tpxyz];
      };
    };
  };


  // fill coord_pmt (pmt->tile number(s))
  Int_t pmt_cnt[25]; 
  for(pp=0; pp<=24; pp++) 
  {
    pmt_cnt[pp]=0; 
    for(int ppp=0; ppp<3; ppp++) coord_pmt[pp][ppp]=0;
  };
  for(tt=1; tt<=36; tt++)
  {
    coord_pmt[GetPMTOfTile(tt)][pmt_cnt[GetPMTOfTile(tt)]] = tt;
    (pmt_cnt[GetPMTOfTile(tt)])++;
  };


  // hexagon side size for event display
  cellsize[kS] = 1;
  cellsize[kL] = 4;

  // line colors
  hexlinecolor[kS] = (Int_t) kGreen+1;
  hexlinecolor[kL] = (Int_t) kAzure+1;


  // initialise event display
  for(sl=0; sl<2; sl++) 
  {
    initTlines = true;
    tile_poly[sl] = this->MakeNewPoly(sl);
    initTlines = false;
    pmt_poly[sl] = this->MakeNewPoly(sl);
    for(tt=0; tt<18; tt++)
    {
      tile_poly[sl]->SetBinContent(coord[sl][tt][kTile]-18*sl,coord[sl][tt][kTile]);
      pmt_poly[sl]->SetBinContent(coord[sl][tt][kTile]-18*sl,coord[sl][tt][kPMT]);
    };
    tile_poly[sl]->SetTitle("BBC tile numbers");
    pmt_poly[sl]->SetTitle("BBC PMT numbers");
  };


  // compute azimuths
  for(tt=1; tt<=36; tt++) azi_tile[tt] = ComputeAzimuthOfTile(tt);
  for(pp=1; pp<=24; pp++) azi_pmt[pp] = ComputeAveAzimuthOfPMT(pp);


  // initialise event polys
  for(sl=0; sl<2; sl++)
  {
    for(ew=0; ew<2; ew++)
    {
      adc_poly[ew][sl] = this->MakeNewPoly(sl);
      tac_poly[ew][sl] = this->MakeNewPoly(sl);
      acc_ev_poly[ew][sl] = this->MakeNewPoly(sl);
    };
  };

  
  // mask PMTs
  for(pp=0; pp<25; pp++) PMTmasked[pp]=false;
  PMTmasked[7]=true;
  PMTmasked[12]=true;
  PMTmasked[18]=true;
  PMTmasked[20]=true;
  PMTmasked[21]=true;
  PMTmasked[22]=true;
  PMTmasked[23]=true;
  PMTmasked[24]=true;


  bbc_canv = new TCanvas("bbc_canv","bbc_canv",1000,500);
  ev_canv = new TCanvas("ev_canv","ev_canv",1000,1000);

  for(sl=0; sl<2; sl++)
  {
    for(ew=0; ew<2; ew++)
    {
      evp_line[ew][sl] = new TLine(0,0,0,0);
      evp_line[ew][sl]->SetLineWidth(2);
      evp_line[ew][sl]->SetLineColor(kMagenta);
      evp_cor_line[ew][sl] = new TLine(0,0,0,0);
      evp_cor_line[ew][sl]->SetLineWidth(2);
      evp_cor_line[ew][sl]->SetLineColor(kOrange-3);
    };
  };

  sigma_mtx = new TMatrix(2,2);
  sigma_evals = new TVector(2);

  ResetEvent(); // zero event variables upon instantiation

  InitRecenteringValues();
};


TH2Poly * BBCtiles::MakeNewPoly(Int_t sl0)
{
  TH2Poly * hc = new TH2Poly();
  for(tt=0; tt<18; tt++)
  {
    xhex = coord[sl0][tt][kX];
    yhex = coord[sl0][tt][kY];
    zhex = coord[sl0][tt][kZ];

    HexToCart(cellsize[sl0],xhex,yhex,zhex,xcart,ycart);

    x[0] = xcart - cellsize[sl0]/2.0;
    x[1] = xcart + cellsize[sl0]/2.0;
    x[2] = xcart + cellsize[sl0];
    x[3] = x[1];
    x[4] = x[0];
    x[5] = xcart - cellsize[sl0];
    
    y[0] = ycart + cellsize[sl0]*sq/2.0;
    y[1] = y[0];
    y[2] = ycart;
    y[3] = ycart - cellsize[sl0]*sq/2.0;
    y[4] = y[3];
    y[5] = y[2];

    for(int kk=0; kk<6; kk++)
    {
      hexline[kk][sl0][tt] = new TLine(x[kk],y[kk],x[(kk+1)%6],y[(kk+1)%6]);
      hexline[kk][sl0][tt]->SetLineColor(hexlinecolor[sl0]);
      hexline[kk][sl0][tt]->SetLineWidth(2);
    };

    hc->AddBin(6, x, y);
  };
  return hc;
};


void BBCtiles::HexToCart(Double_t cellsize0, Int_t x0, Int_t y0, Int_t z0, Double_t &xx0, Double_t &yy0)
{
  xx0 = 1.5 * cellsize0 * x0;
  yy0 = TMath::Sqrt(3)/2.0 * cellsize0 * (y0-z0);
};


/* coordinate accessors */

// return TH2Poly bin number from tile number (to be used with GetSlOfTile)
Int_t BBCtiles::GetBinOfTile(Int_t tile0) { return tile0 - (tile0>18)*18; };

// return sl from tile number
Int_t BBCtiles::GetSlOfTile(Int_t tile0) { return (tile0>18) ? kL:kS; };

// return sl from pmt number
Int_t BBCtiles::GetSlOfPMT(Int_t pmt0) { return (pmt0>16) ? kL:kS; };

// return PMT from tile number
Int_t BBCtiles::GetPMTOfTile(Int_t tile0) { return coord[GetSlOfTile(tile0)][GetBinOfTile(tile0)-1][kPMT]; };

// return tile number(s) from PMT (use whichTile to pick which tile)
Int_t BBCtiles::GetTileOfPMT(Int_t pmt0, Int_t whichTile)
{
  if(whichTile>=0 && whichTile<3)
    return coord_pmt[pmt0][whichTile];
  else return 0;
};

// return tile azimuth
Double_t BBCtiles::GetAzimuthOfTile(Int_t tile0,Bool_t returnDeg)
{
  if(returnDeg==1) return azi_tile[tile0] / TMath::Pi() * 180.0;
  else return azi_tile[tile0];
};
Double_t BBCtiles::GetAveAzimuthOfPMT(Int_t pmt0) { return azi_pmt[pmt0]; };


// return hex coordinates from tile number
Int_t BBCtiles::GetXhexOfTile(Int_t tile0) { return coord[GetSlOfTile(tile0)][GetBinOfTile(tile0)-1][kX]; };
Int_t BBCtiles::GetYhexOfTile(Int_t tile0) { return coord[GetSlOfTile(tile0)][GetBinOfTile(tile0)-1][kY]; };
Int_t BBCtiles::GetZhexOfTile(Int_t tile0) { return coord[GetSlOfTile(tile0)][GetBinOfTile(tile0)-1][kZ]; };


// compute azimuthal angle
Double_t BBCtiles::ComputeAzimuthOfTile(Int_t tile0, Bool_t returnDeg)
{
  Double_t az;
  xhex = GetXhexOfTile(tile0);
  yhex = GetYhexOfTile(tile0);
  zhex = GetZhexOfTile(tile0);
  HexToCart(cellsize[GetSlOfTile(tile0)],xhex,yhex,zhex,xcart,ycart);
  az = TMath::ATan2(ycart,xcart);
  if(returnDeg==1) return az / TMath::Pi() * 180.0; // degrees
  else return az; // radians
};


// compute average azimuthal angle of PMT
Double_t BBCtiles::ComputeAveAzimuthOfPMT(Int_t pmt0)
{
  Int_t tile_current;
  Double_t az_tmp;
  Double_t cos_az=0;
  Double_t sin_az=0;
  Int_t tile_count=0;

  for(int tl=0; tl<3; tl++)
  {
    tile_current = GetTileOfPMT(pmt0,tl);
    if(tile_current>0)
    {
      az_tmp = ComputeAzimuthOfTile(tile_current);
      cos_az += TMath::Cos(az_tmp);
      sin_az += TMath::Sin(az_tmp);
      tile_count++;
    };
  };

  cos_az /= (Double_t)tile_count;
  sin_az /= (Double_t)tile_count;

  return TMath::ATan2(sin_az,cos_az);
};
  

// draw BBC tile and PMT maps
void BBCtiles::DrawBBC()
{
  bbc_canv->Clear();
  bbc_canv->Divide(2,1);
  for(int cc=1; cc<=2; cc++) bbc_canv->GetPad(cc)->SetGrid(1,1);
  bbc_canv->cd(1);
  tile_poly[kL]->Draw("text");
  tile_poly[kS]->Draw("textsame");
  bbc_canv->cd(2);
  pmt_poly[kL]->Draw("text");
  pmt_poly[kS]->Draw("textsame");

  for(int cc=1; cc<=2; cc++)
  {
    bbc_canv->cd(cc);
    for(sl=0; sl<2; sl++)
    {
      for(tt=0; tt<18; tt++)
      {
        for(int xx=0; xx<6; xx++)
        {
          hexline[xx][sl][tt]->Draw();
        };
      };
    };
  };
};


// print BBC tile and PMT maps as well as coordinates
void BBCtiles::PrintBBC()
{
  printf("tile\tPMT\tazimuth(d)\tazimuth(r)\tbin\tsl\ttiles_of_pmt\t<pmt_azimuth>\n");
  printf("----\t---\t----------\t----------\t---\t--\t------------\t-------------\n");
  for(tt=1; tt<=36; tt++)
  {
    printf("%4d\t%3d\t%9.2f\t%9.2f\t%3d\t%2d\t%3d,%3d,%3d\t%13.2f\n",
      tt,
      GetPMTOfTile(tt),
      GetAzimuthOfTile(tt,1),
      GetAzimuthOfTile(tt),
      GetBinOfTile(tt),
      GetSlOfTile(tt),
      GetTileOfPMT(GetPMTOfTile(tt),0),
      GetTileOfPMT(GetPMTOfTile(tt),1),
      GetTileOfPMT(GetPMTOfTile(tt),2),
      GetAveAzimuthOfPMT(GetPMTOfTile(tt))
    );
  };
};


// set new event displays and compute EVPs
void BBCtiles::UpdateEvent()
{
  // reset adc & tac displays
  for(ew=0; ew<2; ew++)
  {
    for(sl=0; sl<2; sl++)
    {
      for(int bb=1; bb<=18; bb++)
      {
        adc_poly[ew][sl]->SetBinContent(bb,0); // TH2poly::Clear doesn't work
        tac_poly[ew][sl]->SetBinContent(bb,0);
      };
      adc_poly[ew][sl]->SetMinimum(0);
      tac_poly[ew][sl]->SetMinimum(0);
      adc_poly[ew][sl]->SetMaximum(4096);
      tac_poly[ew][sl]->SetMaximum(4096);
      EVP[ew][sl]=1000;
      EVP_cor[ew][sl]=1000;
    };
  };

  // update event display
  Int_t tile3,content;
  for(ew=0; ew<2; ew++)
  {
    for(sl=0; sl<2; sl++)
    {
      for(int qq=0; qq<QTN[ew][sl]; qq++)
      {
        for(int t3=0; t3<3; t3++)
        {
          tile3 = GetTileOfPMT(Idx[ew][sl][qq],t3);
          if(tile3>0)
          {
            if(!(PMTmasked[Idx[ew][sl][qq]]))
            {
              content = acc_ev_poly[ew][sl]->GetBinContent(GetBinOfTile(tile3)); // stupid hack to avoid override
              acc_ev_poly[ew][sl]->SetBinContent(GetBinOfTile(tile3),content+1); // " "
              adc_poly[ew][sl]->SetBinContent(GetBinOfTile(tile3), ADC[ew][sl][qq]);
              tac_poly[ew][sl]->SetBinContent(GetBinOfTile(tile3), TAC[ew][sl][qq]);
            };
          };
        };
      };
      ComputeEVP(ew,sl);
      ComputeMoments(ew,sl);
    };
  };
  return;
};


// reset event variables and ADC & TAC displays
void BBCtiles::ResetEvent()
{
  for(sl=0; sl<2; sl++)
  {
    for(ew=0; ew<2; ew++)
    {
      for(int qq=0; qq<16; qq++)
      {
        Idx[ew][sl][qq]=0;
        ADC[ew][sl][qq]=0;
        TAC[ew][sl][qq]=0;
      };
      EVP[ew][sl]=0;
      EVP_cor[ew][sl]=0;
    };
  };
};


// compute event plane (i.e. fourier coefficients)
void BBCtiles::ComputeEVP(Int_t ew0, Int_t sl0)
{
  xfl1 = 0; 
  yfl1 = 0;
  xfl2 = 0; 
  yfl2 = 0;
  Int_t pmt_curr;
  if(QTN[ew0][sl0]>=2)
  {
    for(int qq=0; qq<QTN[ew0][sl0]; qq++)
    {
      pmt_curr = Idx[ew0][sl0][qq];
      if(!(PMTmasked[pmt_curr]))
      {
        xfl1 += ADC[ew0][sl0][qq] * TMath::Cos(GetAveAzimuthOfPMT(pmt_curr));
        yfl1 += ADC[ew0][sl0][qq] * TMath::Sin(GetAveAzimuthOfPMT(pmt_curr));
        xfl2 += ADC[ew0][sl0][qq] * TMath::Cos(2*GetAveAzimuthOfPMT(pmt_curr));
        yfl2 += ADC[ew0][sl0][qq] * TMath::Sin(2*GetAveAzimuthOfPMT(pmt_curr));
      };
    };

    dir[ew0][sl0] = TMath::ATan2(yfl1,xfl1);

    Xflow[ew0][sl0] = xfl2;
    Yflow[ew0][sl0] = yfl2;
    EVP[ew0][sl0] = 0.5 * TMath::ATan2(Yflow[ew0][sl0], Xflow[ew0][sl0]);

    Xflow_cor[ew0][sl0] = 
      ( xfl2 - Xflow_calib[ew0][sl0][QTN[ew0][sl0]-1][0] ) / Xflow_calib[ew0][sl0][QTN[ew0][sl0]-1][1];
    Yflow_cor[ew0][sl0] = 
      ( yfl2 - Yflow_calib[ew0][sl0][QTN[ew0][sl0]-1][0] ) / Yflow_calib[ew0][sl0][QTN[ew0][sl0]-1][1];
    EVP_cor[ew0][sl0] = 0.5 * TMath::ATan2(Yflow_cor[ew0][sl0], Xflow_cor[ew0][sl0]);
  }
  else 
  {
    Xflow[ew0][sl0] = 0;
    Yflow[ew0][sl0] = 0;
    EVP[ew0][sl0] = 1000; // no EVP computed for zero multiplicity

    Xflow_cor[ew0][sl0] = 0;
    Yflow_cor[ew0][sl0] = 0;
    EVP_cor[ew0][sl0] = 1000; // no EVP computed for zero multiplicity
  };
  return;
};


void BBCtiles::ComputeMoments(Int_t ew0, Int_t sl0)
{
  Int_t pmt_curr,tile3,adc_curr;
  Double_t xcart_rot,ycart_rot;
  Double_t xcart_rot_cor,ycart_rot_cor;

  xbar[ew0][sl0]=0;
  ybar[ew0][sl0]=0;
  sigma_x[ew0][sl0]=0;
  sigma_y[ew0][sl0]=0;
  sigma_xy[ew0][sl0]=0;
  sigma_min[ew0][sl0]=0;
  sigma_max[ew0][sl0]=0;
  esum[ew0][sl0]=0;

  left_sum=0;
  right_sum=0;
  directionality[ew0][sl0]=-1;
  left_sum_cor=0;
  right_sum_cor=0;
  directionality_cor[ew0][sl0]=-1;

  if(QTN[ew0][sl0]>0)
  {
    // compute first moments first, which are prerequesite variables to compute second moment
    for(int whichMoment=1; whichMoment<=2; whichMoment++)
    {
      // loop through active channels 
      for(int qq=0; qq<QTN[ew0][sl0]; qq++)
      {
        pmt_curr = Idx[ew0][sl0][qq];
        adc_curr = ADC[ew0][sl0][qq];
        if(!(PMTmasked[pmt_curr]))
        {
          // first obtain tile number; since there is some degeneracy in tile for
          // given PMT, I'm just using the zeroth one for now, since the degenerate
          // channels are masked out for now; no idea what to do otherwise
          tile3 = GetTileOfPMT(pmt_curr,0);

          // get hex coordinates of tile and convert them to cartesian coordinates
          xhex = GetXhexOfTile(tile3);
          yhex = GetYhexOfTile(tile3);
          zhex = GetZhexOfTile(tile3);
          HexToCart(cellsize[GetSlOfTile(tile3)],xhex,yhex,zhex,xcart,ycart);

          // rotate cartesian coordinates by EVP azimuth (for directionality computation)
          xcart_rot =      xcart * TMath::Cos(EVP[ew0][sl0]) + ycart * TMath::Sin(EVP[ew0][sl0]);
          ycart_rot = -1 * xcart * TMath::Sin(EVP[ew0][sl0]) + ycart * TMath::Cos(EVP[ew0][sl0]);
          xcart_rot_cor =      xcart * TMath::Cos(EVP_cor[ew0][sl0]) + ycart * TMath::Sin(EVP_cor[ew0][sl0]);
          ycart_rot_cor = -1 * xcart * TMath::Sin(EVP_cor[ew0][sl0]) + ycart * TMath::Cos(EVP_cor[ew0][sl0]);

          // append to numerators for moments
          if(whichMoment==1)
          {
            xbar[ew0][sl0] += adc_curr * xcart;
            ybar[ew0][sl0] += adc_curr * ycart;

            esum[ew0][sl0] += adc_curr; // only append to esum once

            if(xcart_rot<0) left_sum += adc_curr;
            else if(xcart_rot>0) right_sum += adc_curr;
            if(xcart_rot_cor<0) left_sum_cor += adc_curr;
            else if(xcart_rot_cor>0) right_sum_cor += adc_curr;
          }
          else if(whichMoment==2)
          {
            sigma_x[ew0][sl0]  += adc_curr * (xbar[ew0][sl0]-xcart) * (xbar[ew0][sl0]-xcart);
            sigma_y[ew0][sl0]  += adc_curr * (ybar[ew0][sl0]-ycart) * (ybar[ew0][sl0]-ycart);
            sigma_xy[ew0][sl0] += adc_curr * (xbar[ew0][sl0]-xcart) * (ybar[ew0][sl0]-ycart);
          };
        };
      };

      // finalise moment computations
      if(whichMoment==1)
      {
        xbar[ew0][sl0] /= esum[ew0][sl0];
        ybar[ew0][sl0] /= esum[ew0][sl0];
      }
      else if(whichMoment==2)
      {
        sigma_x[ew0][sl0] /= esum[ew0][sl0];
        sigma_y[ew0][sl0] /= esum[ew0][sl0];
        sigma_xy[ew0][sl0] /= esum[ew0][sl0];
      };
    };

    // compute sigma_min and simga_max
    /* using TMatrix::EigenVectors */
    /*
    (*sigma_mtx)(0,0) = sigma_x[ew0][sl0];
    (*sigma_mtx)(0,1) = sigma_xy[ew0][sl0];
    (*sigma_mtx)(1,0) = sigma_xy[ew0][sl0];
    (*sigma_mtx)(1,1) = sigma_y[ew0][sl0];
    sigma_mtx->EigenVectors(*sigma_evals);
    sigma_min[ew0][sl0] = TMath::Min((*sigma_evals)(0),(*sigma_evals)(1));
    sigma_max[ew0][sl0] = TMath::Max((*sigma_evals)(0),(*sigma_evals)(1));
    */
    /* manually */
    discrim = TMath::Power(sigma_x[ew0][sl0] - sigma_y[ew0][sl0], 2) + 4 * TMath::Power(sigma_xy[ew0][sl0], 2);

    sigma_min_tmp = 0.5 * ( sigma_x[ew0][sl0] + sigma_y[ew0][sl0] - TMath::Sqrt(discrim) );
    sigma_max_tmp = 0.5 * ( sigma_x[ew0][sl0] + sigma_y[ew0][sl0] + TMath::Sqrt(discrim) );
    sigma_min[ew0][sl0] = TMath::Min(sigma_min_tmp,sigma_max_tmp);
    sigma_max[ew0][sl0] = TMath::Max(sigma_min_tmp,sigma_max_tmp);
    
    // compute eigenvector angle.. ambiguity on which one to use
    /* range 0 to +pi */
    //sigma_theta[ew0][sl0] = TMath::ATan2(sigma_max[ew0][sl0]-sigma_x[ew0][sl0], sigma_xy[ew0][sl0]);
    /* range -pi/2 to +pi/2... can impose directionality just like on EVP... */
    sigma_theta[ew0][sl0] = TMath::ATan2(sigma_xy[ew0][sl0], sigma_max[ew0][sl0]-sigma_y[ew0][sl0]);


    // determine directionality 
    if(left_sum == right_sum) directionality[ew0][sl0]=2;
    else if(left_sum < right_sum) directionality[ew0][sl0]=0;
    else directionality[ew0][sl0]=1;

    if(left_sum_cor == right_sum_cor) directionality_cor[ew0][sl0]=2;
    else if(left_sum_cor < right_sum_cor) directionality_cor[ew0][sl0]=0;
    else directionality_cor[ew0][sl0]=1;


    // fix EVP azimuth directionality
    if(directionality[ew0][sl0]==1)
    {
      EVP[ew0][sl0] = (EVP[ew0][sl0]<0) ? EVP[ew0][sl0]+TMath::Pi() : EVP[ew0][sl0]-TMath::Pi();
      sigma_theta[ew0][sl0] = 
        (sigma_theta[ew0][sl0]<0) ? sigma_theta[ew0][sl0]+TMath::Pi() : sigma_theta[ew0][sl0]-TMath::Pi();
    };
    if(directionality_cor[ew0][sl0]==1)
    {
      EVP_cor[ew0][sl0] = (EVP_cor[ew0][sl0]<0) ? EVP_cor[ew0][sl0]+TMath::Pi() : EVP_cor[ew0][sl0]-TMath::Pi();
    };
  };
  return;
};



// draw event
void BBCtiles::DrawEvent()
{
  // EVP lines
  scale = 4*cellsize[kL]; // 1/2 length of evp_line
  for(ew=0; ew<2; ew++)
  {
    for(sl=0; sl<2; sl++)
    {
      if(EVP[ew][sl]<1000)
      {
        switch(directionality[ew][sl])
        {
          case 0:
            evp_line[ew][sl]->SetX1(0);
            evp_line[ew][sl]->SetY1(0);
            evp_line[ew][sl]->SetX2(scale*TMath::Cos(EVP[ew][sl]));
            evp_line[ew][sl]->SetY2(scale*TMath::Sin(EVP[ew][sl]));
            evp_cor_line[ew][sl]->SetX1(0);
            evp_cor_line[ew][sl]->SetY1(0);
            evp_cor_line[ew][sl]->SetX2(scale*TMath::Cos(EVP_cor[ew][sl]));
            evp_cor_line[ew][sl]->SetY2(scale*TMath::Sin(EVP_cor[ew][sl]));
            break;
          case 1:
            evp_line[ew][sl]->SetX1(0);
            evp_line[ew][sl]->SetY1(0);
            evp_line[ew][sl]->SetX2(scale*TMath::Cos(EVP[ew][sl]));
            evp_line[ew][sl]->SetY2(scale*TMath::Sin(EVP[ew][sl]));
            evp_cor_line[ew][sl]->SetX1(0);
            evp_cor_line[ew][sl]->SetY1(0);
            evp_cor_line[ew][sl]->SetX2(scale*TMath::Cos(EVP_cor[ew][sl]));
            evp_cor_line[ew][sl]->SetY2(scale*TMath::Sin(EVP_cor[ew][sl]));
            break;
          case 2:
            evp_line[ew][sl]->SetX1(scale*TMath::Cos(EVP[ew][sl]));
            evp_line[ew][sl]->SetY1(scale*TMath::Sin(EVP[ew][sl]));
            evp_line[ew][sl]->SetX2(-1*scale*TMath::Cos(EVP[ew][sl]));
            evp_line[ew][sl]->SetY2(-1*scale*TMath::Sin(EVP[ew][sl]));
            evp_cor_line[ew][sl]->SetX1(scale*TMath::Cos(EVP_cor[ew][sl]));
            evp_cor_line[ew][sl]->SetY1(scale*TMath::Sin(EVP_cor[ew][sl]));
            evp_cor_line[ew][sl]->SetX2(-1*scale*TMath::Cos(EVP_cor[ew][sl]));
            evp_cor_line[ew][sl]->SetY2(-1*scale*TMath::Sin(EVP_cor[ew][sl]));
            break;
          case -1:
            evp_line[ew][sl]->SetX1(0);
            evp_line[ew][sl]->SetY1(0);
            evp_line[ew][sl]->SetX2(0);
            evp_line[ew][sl]->SetY2(0);
            evp_cor_line[ew][sl]->SetX1(0);
            evp_cor_line[ew][sl]->SetY1(0);
            evp_cor_line[ew][sl]->SetX2(0);
            evp_cor_line[ew][sl]->SetY2(0);
            break;
        };
      }
      else
      {
        evp_line[ew][sl]->SetX1(0);
        evp_line[ew][sl]->SetY1(0);
        evp_line[ew][sl]->SetX2(0);
        evp_line[ew][sl]->SetY2(0);
        evp_cor_line[ew][sl]->SetX1(0);
        evp_cor_line[ew][sl]->SetY1(0);
        evp_cor_line[ew][sl]->SetX2(0);
        evp_cor_line[ew][sl]->SetY2(0);
      }
    };
  };


  // data string
  for(ew=0; ew<2; ew++)
  {
    for(sl=0; sl<2; sl++)
    {
      /*
      datastring[ew][sl] 
        = Form("EVP=%.2f EVP_cor=%.2f #sigma_{x}=%.2f #sigma_{y}=%.2f #sigma_{xy}=%.2f #sigma_{min}=%.2f #sigma_{max}=%.2f",
        EVP[ew][sl],EVP_cor[ew][sl],sigma_x[ew][sl],
        sigma_y[ew][sl],sigma_xy[ew][sl],sigma_min[ew][sl],sigma_max[ew][sl]);
      */
      datastring[ew][sl] 
        = Form("EVP=%.2f EVP_cor=%.2f V%d H%d",
        EVP[ew][sl],EVP_cor[ew][sl],IsVertical(),IsHorizontal());

      adc_poly[ew][sl]->SetTitle(datastring[ew][0].Data()); // for now, no large cells read out, just use small
      tac_poly[ew][sl]->SetTitle("TACs"); // for now, no large cells read out, just use small
    };
  };


  // draw canvas; ADC on top, TAC on bottom; E is left, W is right 
  ev_canv->Clear();
  ev_canv->Divide(2,2);

  for(int ccc=1; ccc<=4; ccc++)
  {
    ev_canv->GetPad(ccc)->SetLogz();
    ev_canv->cd(ccc);
    ew = (ccc-1)%2;
    if(ccc<3)
    {
      adc_poly[ew][kL]->Draw("colz");
      adc_poly[ew][kS]->Draw("colzsame");
    }
    else
    {
      tac_poly[ew][kL]->Draw("colz");
      tac_poly[ew][kS]->Draw("colzsame");
    };

    pmt_poly[kL]->Draw("textsame");
    pmt_poly[kS]->Draw("textsame");

    for(sl=0; sl<2; sl++)
    {
      for(tt=0; tt<18; tt++)
      {
        for(int xx=0; xx<6; xx++)
        {
          hexline[xx][sl][tt]->Draw();
        };
      };
    };

    if(ccc<3) 
    {
      for(sl=0; sl<2; sl++) 
      {
        evp_line[ew][sl]->Draw();
        evp_cor_line[ew][sl]->Draw();
      };
    }
    else 
    {
      for(sl=0; sl<2; sl++) 
      {
        evp_line[ew][sl]->Draw();
        evp_cor_line[ew][sl]->Draw();
      };
    };
  };
};


// EVP class selection
Bool_t BBCtiles::IsVertical()
{
  // vertex cut
  if(TMath::Abs(vertex) > 200) return false;

  div[0] = TMath::Pi()/3.0; // quadrant I bound
  div[1] = 2*TMath::Pi()/3.0; // quadrant II bound
  div[2] = -2*TMath::Pi()/3.0; // quadrant III bound
  div[3] = -1*TMath::Pi()/3.0; // quadrant IV bound

  eee[kE] = EVP_cor[kE][0];
  eee[kW] = EVP_cor[kW][0];

  return ( ( eee[kE] >= div[0] && eee[kE] <= div[1] ) ||
           ( eee[kE] >= div[2] && eee[kE] <= div[3] )
         ) &&
         ( ( eee[kW] >= div[0] && eee[kW] <= div[1] ) ||
           ( eee[kW] >= div[2] && eee[kW] <= div[3] )
         );
};

Bool_t BBCtiles::IsHorizontal()
{
  // vertex cut
  if(TMath::Abs(vertex) > 200) return false;

  div[0] = TMath::Pi()/6.0; // quadrant I bound
  div[1] = 5*TMath::Pi()/6.0; // quadrant II bound
  div[2] = -5*TMath::Pi()/6.0; // quadrant III bound
  div[3] = -1*TMath::Pi()/6.0; // quadrant IV bound

  eee[kE] = EVP_cor[kE][0];
  eee[kW] = EVP_cor[kW][0];

  return ( ( eee[kE] >= div[3] && eee[kE] <= div[0] ) ||
           ( eee[kE] >= div[1] && eee[kE] <= TMath::Pi() ) ||
           ( eee[kE] >= -1*TMath::Pi() && eee[kE] <= div[2] )
         ) &&
         ( ( eee[kW] >= div[3] && eee[kW] <= div[0] ) ||
           ( eee[kW] >= div[1] && eee[kW] <= TMath::Pi() ) ||
           ( eee[kW] >= -1*TMath::Pi() && eee[kW] <= div[2] )
         );
};


// initialize flow vector means & rms for re-centering correction
// (copy and paste output from EVPdiagnostics.C)
void BBCtiles::InitRecenteringValues()
{
  Xflow_calib[0][0][0][0]=46.30; Xflow_calib[0][0][0][1]=383.81;
  Yflow_calib[0][0][0][0]=60.47; Yflow_calib[0][0][0][1]=272.71;
  Xflow_calib[0][0][1][0]=41.80; Xflow_calib[0][0][1][1]=497.64;
  Yflow_calib[0][0][1][0]=68.18; Yflow_calib[0][0][1][1]=348.97;
  Xflow_calib[0][0][2][0]=33.54; Xflow_calib[0][0][2][1]=597.32;
  Yflow_calib[0][0][2][0]=66.66; Yflow_calib[0][0][2][1]=426.96;
  Xflow_calib[0][0][3][0]=14.62; Xflow_calib[0][0][3][1]=682.20;
  Yflow_calib[0][0][3][0]=54.61; Yflow_calib[0][0][3][1]=483.60;
  Xflow_calib[0][0][4][0]=-5.41; Xflow_calib[0][0][4][1]=768.23;
  Yflow_calib[0][0][4][0]=38.44; Yflow_calib[0][0][4][1]=558.66;
  Xflow_calib[0][0][5][0]=-18.05; Xflow_calib[0][0][5][1]=849.70;
  Yflow_calib[0][0][5][0]=26.01; Yflow_calib[0][0][5][1]=620.26;
  Xflow_calib[0][0][6][0]=-29.69; Xflow_calib[0][0][6][1]=922.55;
  Yflow_calib[0][0][6][0]=-1.77; Yflow_calib[0][0][6][1]=680.60;
  Xflow_calib[0][0][7][0]=-61.19; Xflow_calib[0][0][7][1]=1006.72;
  Yflow_calib[0][0][7][0]=-33.90; Yflow_calib[0][0][7][1]=753.23;
  Xflow_calib[0][0][8][0]=-80.92; Xflow_calib[0][0][8][1]=1087.05;
  Yflow_calib[0][0][8][0]=-84.85; Yflow_calib[0][0][8][1]=818.35;
  Xflow_calib[0][0][9][0]=-110.94; Xflow_calib[0][0][9][1]=1180.75;
  Yflow_calib[0][0][9][0]=-133.45; Yflow_calib[0][0][9][1]=890.75;
  Xflow_calib[0][0][10][0]=-131.77; Xflow_calib[0][0][10][1]=1274.49;
  Yflow_calib[0][0][10][0]=-197.75; Yflow_calib[0][0][10][1]=988.28;
  Xflow_calib[0][0][11][0]=-172.23; Xflow_calib[0][0][11][1]=1401.41;
  Yflow_calib[0][0][11][0]=-268.87; Yflow_calib[0][0][11][1]=1081.05;
  Xflow_calib[0][0][12][0]=-182.98; Xflow_calib[0][0][12][1]=1523.05;
  Yflow_calib[0][0][12][0]=-355.78; Yflow_calib[0][0][12][1]=1188.18;
  Xflow_calib[0][0][13][0]=-191.81; Xflow_calib[0][0][13][1]=1679.05;
  Yflow_calib[0][0][13][0]=-434.86; Yflow_calib[0][0][13][1]=1316.19;
  Xflow_calib[0][0][14][0]=-216.74; Xflow_calib[0][0][14][1]=1833.76;
  Yflow_calib[0][0][14][0]=-507.84; Yflow_calib[0][0][14][1]=1453.85;
  Xflow_calib[0][0][15][0]=-192.05; Xflow_calib[0][0][15][1]=1966.15;
  Yflow_calib[0][0][15][0]=-600.54; Yflow_calib[0][0][15][1]=1533.06;
  Xflow_calib[0][1][0][0]=0.00; Xflow_calib[0][1][0][1]=1.00;
  Yflow_calib[0][1][0][0]=0.00; Yflow_calib[0][1][0][1]=1.00;
  Xflow_calib[0][1][1][0]=0.00; Xflow_calib[0][1][1][1]=1.00;
  Yflow_calib[0][1][1][0]=0.00; Yflow_calib[0][1][1][1]=1.00;
  Xflow_calib[0][1][2][0]=0.00; Xflow_calib[0][1][2][1]=1.00;
  Yflow_calib[0][1][2][0]=0.00; Yflow_calib[0][1][2][1]=1.00;
  Xflow_calib[0][1][3][0]=0.00; Xflow_calib[0][1][3][1]=1.00;
  Yflow_calib[0][1][3][0]=0.00; Yflow_calib[0][1][3][1]=1.00;
  Xflow_calib[0][1][4][0]=0.00; Xflow_calib[0][1][4][1]=1.00;
  Yflow_calib[0][1][4][0]=0.00; Yflow_calib[0][1][4][1]=1.00;
  Xflow_calib[0][1][5][0]=0.00; Xflow_calib[0][1][5][1]=1.00;
  Yflow_calib[0][1][5][0]=0.00; Yflow_calib[0][1][5][1]=1.00;
  Xflow_calib[0][1][6][0]=0.00; Xflow_calib[0][1][6][1]=1.00;
  Yflow_calib[0][1][6][0]=0.00; Yflow_calib[0][1][6][1]=1.00;
  Xflow_calib[0][1][7][0]=0.00; Xflow_calib[0][1][7][1]=1.00;
  Yflow_calib[0][1][7][0]=0.00; Yflow_calib[0][1][7][1]=1.00;
  Xflow_calib[0][1][8][0]=0.00; Xflow_calib[0][1][8][1]=1.00;
  Yflow_calib[0][1][8][0]=0.00; Yflow_calib[0][1][8][1]=1.00;
  Xflow_calib[0][1][9][0]=0.00; Xflow_calib[0][1][9][1]=1.00;
  Yflow_calib[0][1][9][0]=0.00; Yflow_calib[0][1][9][1]=1.00;
  Xflow_calib[0][1][10][0]=0.00; Xflow_calib[0][1][10][1]=1.00;
  Yflow_calib[0][1][10][0]=0.00; Yflow_calib[0][1][10][1]=1.00;
  Xflow_calib[0][1][11][0]=0.00; Xflow_calib[0][1][11][1]=1.00;
  Yflow_calib[0][1][11][0]=0.00; Yflow_calib[0][1][11][1]=1.00;
  Xflow_calib[0][1][12][0]=0.00; Xflow_calib[0][1][12][1]=1.00;
  Yflow_calib[0][1][12][0]=0.00; Yflow_calib[0][1][12][1]=1.00;
  Xflow_calib[0][1][13][0]=0.00; Xflow_calib[0][1][13][1]=1.00;
  Yflow_calib[0][1][13][0]=0.00; Yflow_calib[0][1][13][1]=1.00;
  Xflow_calib[0][1][14][0]=0.00; Xflow_calib[0][1][14][1]=1.00;
  Yflow_calib[0][1][14][0]=0.00; Yflow_calib[0][1][14][1]=1.00;
  Xflow_calib[0][1][15][0]=0.00; Xflow_calib[0][1][15][1]=1.00;
  Yflow_calib[0][1][15][0]=0.00; Yflow_calib[0][1][15][1]=1.00;
  Xflow_calib[1][0][0][0]=77.72; Xflow_calib[1][0][0][1]=314.81;
  Yflow_calib[1][0][0][0]=35.01; Yflow_calib[1][0][0][1]=318.92;
  Xflow_calib[1][0][1][0]=125.10; Xflow_calib[1][0][1][1]=426.76;
  Yflow_calib[1][0][1][0]=65.45; Yflow_calib[1][0][1][1]=388.24;
  Xflow_calib[1][0][2][0]=183.30; Xflow_calib[1][0][2][1]=531.37;
  Yflow_calib[1][0][2][0]=96.34; Yflow_calib[1][0][2][1]=473.36;
  Xflow_calib[1][0][3][0]=226.75; Xflow_calib[1][0][3][1]=628.88;
  Yflow_calib[1][0][3][0]=137.14; Yflow_calib[1][0][3][1]=538.70;
  Xflow_calib[1][0][4][0]=262.33; Xflow_calib[1][0][4][1]=702.26;
  Yflow_calib[1][0][4][0]=165.09; Yflow_calib[1][0][4][1]=594.98;
  Xflow_calib[1][0][5][0]=278.42; Xflow_calib[1][0][5][1]=771.01;
  Yflow_calib[1][0][5][0]=188.92; Yflow_calib[1][0][5][1]=649.69;
  Xflow_calib[1][0][6][0]=293.18; Xflow_calib[1][0][6][1]=850.79;
  Yflow_calib[1][0][6][0]=203.07; Yflow_calib[1][0][6][1]=715.43;
  Xflow_calib[1][0][7][0]=275.24; Xflow_calib[1][0][7][1]=919.95;
  Yflow_calib[1][0][7][0]=209.23; Yflow_calib[1][0][7][1]=776.18;
  Xflow_calib[1][0][8][0]=261.46; Xflow_calib[1][0][8][1]=994.65;
  Yflow_calib[1][0][8][0]=205.10; Yflow_calib[1][0][8][1]=834.35;
  Xflow_calib[1][0][9][0]=252.59; Xflow_calib[1][0][9][1]=1058.39;
  Yflow_calib[1][0][9][0]=181.52; Yflow_calib[1][0][9][1]=905.22;
  Xflow_calib[1][0][10][0]=225.24; Xflow_calib[1][0][10][1]=1131.57;
  Yflow_calib[1][0][10][0]=182.74; Yflow_calib[1][0][10][1]=966.80;
  Xflow_calib[1][0][11][0]=179.09; Xflow_calib[1][0][11][1]=1245.70;
  Yflow_calib[1][0][11][0]=136.07; Yflow_calib[1][0][11][1]=1045.91;
  Xflow_calib[1][0][12][0]=151.53; Xflow_calib[1][0][12][1]=1370.14;
  Yflow_calib[1][0][12][0]=114.35; Yflow_calib[1][0][12][1]=1149.54;
  Xflow_calib[1][0][13][0]=80.88; Xflow_calib[1][0][13][1]=1501.78;
  Yflow_calib[1][0][13][0]=56.89; Yflow_calib[1][0][13][1]=1250.89;
  Xflow_calib[1][0][14][0]=64.02; Xflow_calib[1][0][14][1]=1628.06;
  Yflow_calib[1][0][14][0]=-3.22; Yflow_calib[1][0][14][1]=1353.10;
  Xflow_calib[1][0][15][0]=85.98; Xflow_calib[1][0][15][1]=1768.64;
  Yflow_calib[1][0][15][0]=-112.77; Yflow_calib[1][0][15][1]=1465.44;
  Xflow_calib[1][1][0][0]=0.00; Xflow_calib[1][1][0][1]=1.00;
  Yflow_calib[1][1][0][0]=0.00; Yflow_calib[1][1][0][1]=1.00;
  Xflow_calib[1][1][1][0]=0.00; Xflow_calib[1][1][1][1]=1.00;
  Yflow_calib[1][1][1][0]=0.00; Yflow_calib[1][1][1][1]=1.00;
  Xflow_calib[1][1][2][0]=0.00; Xflow_calib[1][1][2][1]=1.00;
  Yflow_calib[1][1][2][0]=0.00; Yflow_calib[1][1][2][1]=1.00;
  Xflow_calib[1][1][3][0]=0.00; Xflow_calib[1][1][3][1]=1.00;
  Yflow_calib[1][1][3][0]=0.00; Yflow_calib[1][1][3][1]=1.00;
  Xflow_calib[1][1][4][0]=0.00; Xflow_calib[1][1][4][1]=1.00;
  Yflow_calib[1][1][4][0]=0.00; Yflow_calib[1][1][4][1]=1.00;
  Xflow_calib[1][1][5][0]=0.00; Xflow_calib[1][1][5][1]=1.00;
  Yflow_calib[1][1][5][0]=0.00; Yflow_calib[1][1][5][1]=1.00;
  Xflow_calib[1][1][6][0]=0.00; Xflow_calib[1][1][6][1]=1.00;
  Yflow_calib[1][1][6][0]=0.00; Yflow_calib[1][1][6][1]=1.00;
  Xflow_calib[1][1][7][0]=0.00; Xflow_calib[1][1][7][1]=1.00;
  Yflow_calib[1][1][7][0]=0.00; Yflow_calib[1][1][7][1]=1.00;
  Xflow_calib[1][1][8][0]=0.00; Xflow_calib[1][1][8][1]=1.00;
  Yflow_calib[1][1][8][0]=0.00; Yflow_calib[1][1][8][1]=1.00;
  Xflow_calib[1][1][9][0]=0.00; Xflow_calib[1][1][9][1]=1.00;
  Yflow_calib[1][1][9][0]=0.00; Yflow_calib[1][1][9][1]=1.00;
  Xflow_calib[1][1][10][0]=0.00; Xflow_calib[1][1][10][1]=1.00;
  Yflow_calib[1][1][10][0]=0.00; Yflow_calib[1][1][10][1]=1.00;
  Xflow_calib[1][1][11][0]=0.00; Xflow_calib[1][1][11][1]=1.00;
  Yflow_calib[1][1][11][0]=0.00; Yflow_calib[1][1][11][1]=1.00;
  Xflow_calib[1][1][12][0]=0.00; Xflow_calib[1][1][12][1]=1.00;
  Yflow_calib[1][1][12][0]=0.00; Yflow_calib[1][1][12][1]=1.00;
  Xflow_calib[1][1][13][0]=0.00; Xflow_calib[1][1][13][1]=1.00;
  Yflow_calib[1][1][13][0]=0.00; Yflow_calib[1][1][13][1]=1.00;
  Xflow_calib[1][1][14][0]=0.00; Xflow_calib[1][1][14][1]=1.00;
  Yflow_calib[1][1][14][0]=0.00; Yflow_calib[1][1][14][1]=1.00;
  Xflow_calib[1][1][15][0]=0.00; Xflow_calib[1][1][15][1]=1.00;
  Yflow_calib[1][1][15][0]=0.00; Yflow_calib[1][1][15][1]=1.00;
};
