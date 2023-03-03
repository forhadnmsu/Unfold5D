#include <iostream>
#include <sstream>
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include <TStyle.h>
#include"TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include <ROOT/RDataFrame.hxx>
#include "ana_const.h"
using namespace std;
void MakeTestHist5D(){
	gStyle->SetOptFit(1);
	TH2::SetDefaultSumw2();
	ROOT::RDataFrame dmyDf_mc1_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test1.root"); 
	ROOT::RDataFrame dmyDf_mc2_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test2.root");
	ROOT::RDataFrame dmyDf_mc3_("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Test3.root");

	auto ReWeight = [](float cosThetaCS, float PhiCS, float wt) { return wt*(1+ 0.9*cosThetaCS*cosThetaCS+(2*0.2)*(sqrt(1-cosThetaCS*cosThetaCS)*cosThetaCS*cos(PhiCS))+0.2*(1-cosThetaCS*cosThetaCS)*cos(2*PhiCS)*0.5)/(1+ cosThetaCS*cosThetaCS);};

	auto dmyDf_mc1 = dmyDf_mc1_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});
	auto dmyDf_mc2 = dmyDf_mc2_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});
	auto dmyDf_mc3 = dmyDf_mc3_.Define("ParReweight", ReWeight, {"true_costh","true_phi","weight"});

	int events= 50000;
	int steps1=*dmyDf_mc1.Count()/events;
	int steps2=*dmyDf_mc2.Count()/events;
	int steps3=*dmyDf_mc3.Count()/events;
	int range1 = *dmyDf_mc1.Count()/steps1;
	int range2 = *dmyDf_mc2.Count()/steps2;
	int range3 = *dmyDf_mc3.Count()/steps3;

	TFile *myfile = new TFile("testMC2D.root","recreate");
	int split=0;

	TH1D* h_5D_xF = new TH1D("h_5D_xF", "",nBinxF,xFMin,xFMax);
	TH1D* h_5D_pT = new TH1D("h_5D_pT", "",nBinpT, pTMin, pTMax);
	TH1D* h_5D_mass = new TH1D("h_5D_mass", "",nBinMass, MassMin, MassMax);
	TH1D* h_5D_costh = new TH1D("h_5D_costh", "",nBinCosth, costhMin, costhMax);
	TH1D* h_5D_phi = new TH1D("h_5D_phi", "",nBinPhi, phiMin, phiMax);

	int nbinX = h_5D_xF->GetNbinsX();
	int nbinY = h_5D_pT->GetNbinsX();
	int nbinZ = h_5D_mass->GetNbinsX();
	int nbinT = h_5D_phi->GetNbinsX();
	int nbinU = h_5D_costh->GetNbinsX();
	int nbinx = nbinX*nbinY*nbinZ*nbinT*nbinU;
	TH1D* h_X_5D[100];

	 for (int i_split = 1; i_split <=steps1 ; i_split++){
		 int begin_= 1+(i_split-1)*range1;
		 int end_= range1* i_split;
		 h_X_5D[i_split-1] = new TH1D(Form("h_X_5D_%i",i_split), "",nbinx,1, nbinx);

		 dmyDf_mc1.Range(begin_,end_).Foreach([i_split, h_X_5D, h_5D_xF, h_5D_pT, h_5D_mass, h_5D_costh, h_5D_phi](float reco_mass, float tr_mass, float reco_phi, float tr_phi, float reco_costh, float tr_costh, float reco_xF,  float tr_xF, float reco_pT, float tr_pT, float wt,double inj_wt, int trig){  
				 cout << "file number in root: "<< i_split << endl;
				 int x_mass = h_5D_mass->GetXaxis()->FindBin(tr_mass);
				 cout << "mass:" << tr_mass << "bin: "<< x_mass << endl;
				 int nbinX = h_5D_xF->GetNbinsX();
				 int nbinY = h_5D_pT->GetNbinsX();
				 int nbinZ = h_5D_mass->GetNbinsX();
				 int nbinT = h_5D_phi->GetNbinsX();
				 int nbinU = h_5D_costh->GetNbinsX();
				 int i_xF_= h_5D_xF->GetXaxis()->FindBin(reco_xF);
				 int i_pT_= h_5D_pT->GetXaxis()->FindBin(reco_pT);
				 int i_mass_= h_5D_mass->GetXaxis()->FindBin(reco_mass);
				 int i_phi_= h_5D_phi->GetXaxis()->FindBin(reco_phi);
				 int i_costh_= h_5D_costh->GetXaxis()->FindBin(reco_costh);
				 int n1 = h_5D_xF->GetNbinsX();
				 int n2 = h_5D_pT->GetNbinsX();
				 int n3 = h_5D_mass->GetNbinsX();
				 int n4 = h_5D_phi->GetNbinsX();
				 int n5 = h_5D_costh->GetNbinsX();
				 int BinIndex1D_reco = (i_xF_-1) * (n2 * n3 * n4 * n5) + (i_pT_-1) * (n3 * n4 * n5) + (i_mass_-1) * (n4 * n5) + (i_phi_-1) *n5 + (i_costh_-1)+1;
				 //h_X_5D[i_split-1]->Fill(BinIndex1D_reco, wt);
				 h_X_5D[i_split-1]->Fill(BinIndex1D_reco, inj_wt);
				 }, {"mass","true_mass", "phi", "true_phi", "costh", "true_costh", "xF", "true_xF", "pT", "true_pT", "weight", "ParReweight","fpga1"} );

		 char name[100]; 
		 sprintf(name, "%s%i","costh_test_",i_split);
		 TH2D* h = (TH2D*)h_X_5D[i_split-1]->Clone();    
		 h->SetName(name); h->Write(name,TObject::kWriteDelete);
	}

	//for file 2
	  for (int i_split = 1; i_split <=steps2 ; i_split++){
                 int begin_= 1+(i_split-1)*range2;
                 int end_= range2* i_split;
                 h_X_5D[i_split-1+steps1] = new TH1D(Form("h_X_5D_%i",i_split+steps1), "",nbinx,1, nbinx);

                 dmyDf_mc2.Range(begin_,end_).Foreach([i_split,steps1, h_X_5D, h_5D_xF, h_5D_pT, h_5D_mass, h_5D_costh, h_5D_phi](float reco_mass, float tr_mass, float reco_phi, float tr_phi, float reco_costh, float tr_costh, float reco_xF,  float tr_xF, float reco_pT, float tr_pT, float wt,double inj_wt, int trig){  
				 cout << "file number in root: "<< i_split+steps1 << endl;
                                 int x_mass = h_5D_mass->GetXaxis()->FindBin(tr_mass);
                                 cout << "mass:" << tr_mass << "bin: "<< x_mass << endl;
                                 int nbinX = h_5D_xF->GetNbinsX();
                                 int nbinY = h_5D_pT->GetNbinsX();
                                 int nbinZ = h_5D_mass->GetNbinsX();
                                 int nbinT = h_5D_phi->GetNbinsX();
                                 int nbinU = h_5D_costh->GetNbinsX();
                                 int i_xF_= h_5D_xF->GetXaxis()->FindBin(reco_xF);
                                 int i_pT_= h_5D_pT->GetXaxis()->FindBin(reco_pT);
                                 int i_mass_= h_5D_mass->GetXaxis()->FindBin(reco_mass);
                                 int i_phi_= h_5D_phi->GetXaxis()->FindBin(reco_phi);
                                 int i_costh_= h_5D_costh->GetXaxis()->FindBin(reco_costh);
                                 int n1 = h_5D_xF->GetNbinsX();
                                 int n2 = h_5D_pT->GetNbinsX();
                                 int n3 = h_5D_mass->GetNbinsX();
                                 int n4 = h_5D_phi->GetNbinsX();
                                 int n5 = h_5D_costh->GetNbinsX();
                                 int BinIndex1D_reco = (i_xF_-1) * (n2 * n3 * n4 * n5) + (i_pT_-1) * (n3 * n4 * n5) + (i_mass_-1) * (n4 * n5) + (i_phi_-1) *n5 + (i_costh_-1)+1;
                                 //h_X_5D[i_split-1+steps1]->Fill(BinIndex1D_reco, wt);
                                 h_X_5D[i_split-1+steps1]->Fill(BinIndex1D_reco, inj_wt);
				 }, {"mass","true_mass", "phi", "true_phi", "costh", "true_costh", "xF", "true_xF", "pT", "true_pT", "weight", "ParReweight","fpga1"} );

                 char name[100]; 
                 sprintf(name, "%s%i","costh_test_",i_split+steps1);
                 TH2D* h = (TH2D*)h_X_5D[i_split-1+steps1]->Clone();    
                 h->SetName(name); h->Write(name,TObject::kWriteDelete);
        } 

	//for file 3
          for (int i_split = 1; i_split <=steps3 ; i_split++){
                 int begin_= 1+(i_split-1)*range3;
                 int end_= range3* i_split;
                 h_X_5D[i_split-1+steps1+steps2] = new TH1D(Form("h_X_5D_%i",i_split+steps1+steps2), "",nbinx,1, nbinx);

                 dmyDf_mc3.Range(begin_,end_).Foreach([i_split,steps1, steps2, h_X_5D, h_5D_xF, h_5D_pT, h_5D_mass, h_5D_costh, h_5D_phi](float reco_mass, float tr_mass, float reco_phi, float tr_phi, float reco_costh, float tr_costh, float reco_xF,  float tr_xF, float reco_pT, float tr_pT, float wt, double inj_wt,int trig){  
				 cout << "file number in root: "<< i_split+steps1+steps2 << endl;
                                 int x_mass = h_5D_mass->GetXaxis()->FindBin(tr_mass);
                                 cout << "mass:" << tr_mass << "bin: "<< x_mass << endl;
                                 int nbinX = h_5D_xF->GetNbinsX();
                                 int nbinY = h_5D_pT->GetNbinsX();
                                 int nbinZ = h_5D_mass->GetNbinsX();
                                 int nbinT = h_5D_phi->GetNbinsX();
                                 int nbinU = h_5D_costh->GetNbinsX();
                                 int i_xF_= h_5D_xF->GetXaxis()->FindBin(reco_xF);
                                 int i_pT_= h_5D_pT->GetXaxis()->FindBin(reco_pT);
                                 int i_mass_= h_5D_mass->GetXaxis()->FindBin(reco_mass);
                                 int i_phi_= h_5D_phi->GetXaxis()->FindBin(reco_phi);
                                 int i_costh_= h_5D_costh->GetXaxis()->FindBin(reco_costh);
                                 int n1 = h_5D_xF->GetNbinsX();
                                 int n2 = h_5D_pT->GetNbinsX();
                                 int n3 = h_5D_mass->GetNbinsX();
                                 int n4 = h_5D_phi->GetNbinsX();
                                 int n5 = h_5D_costh->GetNbinsX();
                                 int BinIndex1D_reco = (i_xF_-1) * (n2 * n3 * n4 * n5) + (i_pT_-1) * (n3 * n4 * n5) + (i_mass_-1) * (n4 * n5) + (i_phi_-1) *n5 + (i_costh_-1)+1;
                                 //h_X_5D[i_split-1+steps1+steps2]->Fill(BinIndex1D_reco, wt);
                                 h_X_5D[i_split-1+steps1+steps2]->Fill(BinIndex1D_reco, inj_wt);
				 }, {"mass","true_mass", "phi", "true_phi", "costh", "true_costh", "xF", "true_xF", "pT", "true_pT", "weight", "ParReweight","fpga1"} );

                 char name[100]; 
                 sprintf(name, "%s%i","costh_test_",i_split+steps1+steps2);
                 TH2D* h = (TH2D*)h_X_5D[i_split-1+steps1+steps2]->Clone();    
                 h->SetName(name); h->Write(name,TObject::kWriteDelete);
        }   

//end files

	myfile->Close();
	myfile->ls();
}
