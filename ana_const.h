#ifndef _ANA_CONSTANTS_H
#define _ANA_CONSTANTS_H	
//global variables and binning

const int nBinMass =2; 
const double MassMin =5.0;
const double MassMax =7.0;

const int nBinxF =2;  
const double xFMin =0.0;  
const double xFMax =1.0; 

const int nBinpT =1; 
const double pTMin =0.0; 
const double pTMax =2.0; 

const int nBinCosth =16;
const double costhMin= -0.9;
const double costhMax= 0.9;
const int nBinPhi = 10; 
const double phiMin = -M_PI;
const double phiMax = M_PI;
const double cos_low=-0.4;
const double cos_high= 0.4;

string posCut = "atan(y1_st3/y1_st1)>1.086 && atan(y1_st3/y1_st1) <1.562&& z1_v > -45.0  && z1_v < 125.0&& x1_t*x1_t+(y1_t-1.6)*(y1_t-1.6) <1425.0 && x1_d*x1_d+(y1_d -1.6)*(y1_d-1.6) <141.0&& abs(px1_st1-px1_st3)> 0.413 &&  abs(px1_st1-px1_st3)<0.420 && abs(py1_st1-py1_st3)< 0.005&& abs(pz1_st1-pz1_st3)<0.074&&(chisq1_target-chisq1_dump) >-2.5 && (chisq1_target-chisq1_dump)<77.5&& pz1_st1 > 10 && pz1_st1 < 78 && nHits1 > 13 &&(y1_st1*y1_st3) > 0 && chisq1_dump/chisq1_upstream < 0.165&&";

string negCut ="atan(y2_st3/y2_st1)>0.666 && atan(y2_st3/y2_st1) <1.506&& z2_v > -45.0 && z2_v < 135.0&& x2_t*x2_t+(y2_t-1.6)*(y2_t-1.6) <1485.0 && x2_d*x2_d+(y2_d-1.6)*(y2_d-1.6) <201.0&& abs(px2_st1-px2_st3)> 0.414 &&  abs(px2_st1-px2_st3)<0.420 && abs(py2_st1-py2_st3)< 0.005&& abs(pz2_st1-pz2_st3)<0.05&&(chisq2_target-chisq2_dump) >-2.5 && (chisq2_target-chisq2_dump)<77.5&&  pz2_st1 > 9 && pz2_st1 < 78 && nHits2 > 13 &&(y2_st1*y2_st3) > 0&& chisq2_dump/chisq2_upstream < 0.195 &&";
 string dimuCut= "dz>0.0 && dz< 1190.0&& abs(dx) < 0.262 && abs(dy-1.6) < 0.225 &&  abs(dpx) < 1.675 && abs(dpy) < 2.4 && dpz > 40.0 && dpz < 114.0 && abs(x1_st1 + x2_st1) < 30 && dpx*dpx + dpy*dpy < 5.3 && dx*dx + (dy-1.6)*(dy-1.6) < 0.105  && y1_st3*y2_st3 < 0 && nHits1 + nHits2 > 29 && nHits1St1 + nHits2St1 >9 &&chisq_dimuon>0.0&&chisq_dimuon<11.5&&(chisq1_dump+chisq2_dump-chisq_dimuon)<19.27&&abs(z1_v - z2_v)<125.0&&";

string physics = "xF > -0.18 && xF <0.94&& xT > 0.007 && xT <0.54 && mass>5. && mass <8.0 && pT <2.0 && D1<200";
string trigger = "fpga1==1&&";
//============Final Cut
string cutTrue = "true_pT<2.0&&true_mass>5. && true_mass<7.0 &&true_xF>0.0 && true_xF<1.0";
string cutRecoData = posCut+negCut+dimuCut+physics; 
string cutRecoMC= cutTrue+ "&&"+trigger +cutRecoData;

//string cutRecoMC = cutTrue+"&&pT<2.0&&mass>5. && mass<7.0&&fpga1==1&&D1<200 &&xF>0.0 && xF<1.";

TH1D * Make2Dto1D(
                        const char* hname,
                         TH2D* h_phi_costh)
{
        int nbinX = h_phi_costh->GetNbinsX();
        int nbinY = h_phi_costh->GetNbinsY();
        int nbin = nbinX*nbinY;
        //TH1D* h_phi_costh_1D = new TH1D("h_phi_costh_1D","h_phi_costh_1D",nbin, 1+0.5,nbin+0.5);
        TH1D* h_phi_costh_1D = new TH1D("h_phi_costh_1D","h_phi_costh_1D",nbin, 1,nbin);
        for (int ii=1; ii<= nbinX; ii++){
                for (int jj=1; jj<= nbinY; jj++){
                        //cout << "get bin number from TH2: \t" << "ii: "<< ii << " jj "<< jj<< " entry bin: "<< h_phi_costh->GetBin(ii,jj) << "we set: "<< ((ii-1)*nbinY + jj)<<endl;

                        int BinIndex1D =  ((ii-1)*nbinY + jj);  //id bin number
                        float x_2d = h_phi_costh->GetBinContent(ii,jj);
                        float x_2d_error = h_phi_costh->GetBinError(ii,jj);
                        h_phi_costh_1D->SetBinContent(BinIndex1D,x_2d);
                        h_phi_costh_1D->SetBinError(BinIndex1D,x_2d_error);
                }
        }
        return h_phi_costh_1D;
}

TH2D * Make1Dto2D(
                        const char* hname,
                        TH1D* hX, TH2D* phi_costh)
{
        TH2D*h = (TH2D*)phi_costh->Clone();
        int nbinX = phi_costh->GetNbinsX();
        int nbinY = phi_costh->GetNbinsY();
        for (int ii=1; ii<= nbinX; ii++){
                for (int jj=1; jj<= nbinY; jj++){

                        int index=  (ii-1)*nbinY + jj;
                        h->SetBinContent(ii,jj,hX->GetBinContent(index));
                        h->SetBinError(ii,jj,hX->GetBinError(index));
                }
        }
        return h;
}

TH2D * ResetEmptyBins(
                        const char* hname,
                         TH2D* h_phi_costh )
 {
        int nbinX = h_phi_costh->GetNbinsX();
        int nbinY = h_phi_costh->GetNbinsY();
        TH2D *h = (TH2D*)h_phi_costh->Clone();
        for (int ii=1; ii<= nbinX; ii++){
                for (int jj=1; jj<= nbinY; jj++){

                        if(h_phi_costh->GetBinContent(ii,jj) <2){
                                h->SetBinContent(ii,jj,0);
                                h->SetBinError(ii,jj,0);
                        };
                }
        }
        return h;
}
#endif
