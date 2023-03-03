#include "ana_const.h"
#include "histogram_utils.h"
R__LOAD_LIBRARY(/Users/forhadhossain/test/RooUnfold/build/libRooUnfold.dylib)
void ana(){
	gStyle->SetOptStat(1111111);
	gStyle->SetStatY(0.9);                
	gStyle->SetStatX(0.9);                
	gStyle->SetStatW(0.4);                
	gStyle->SetStatH(0.2);                

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	//gStyle->SetOptStat(0);
	ROOT::RDataFrame Df_mc1("result_mc", "data/DY_4pi_GMC_Jan03_LD2_Train.root"); 
	string passed = cutRecoMC;
	auto Df_train = Df_mc1.Define("passed", passed);

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
	RooUnfoldResponse response(nbinx,1, nbinx);
	int bin1Dx= nbinT*nbinU;
	RooUnfoldResponse response1D(bin1Dx,1, bin1Dx);
	TH2D* phi_theta_data_tr = new TH2D("phi_theta_data_tr", "Data True  #phi_{CS} Vs. cos#theta_{CS}",nBinPhi,    phiMin, phiMax, nBinCosth, costhMin, costhMax);                 
	TH2D* phi_theta_gen = new TH2D("phi_theta_gen", "MC Generated #phi_{CS} Vs. cos#theta_{CS};#phi_{CS};cos#theta_{CS}",nBinPhi,phiMin,phiMax,nBinCosth,costhMin,costhMax);
	TH2D* phi_theta_reco = new TH2D("phi_theta_reco", "MC Reconstructed  #phi_{CS} Vs. cos#theta_{CS};#phi_{CS};cos#theta_{CS}", nBinPhi,phiMin,phiMax,nBinCosth,costhMin,costhMax);
	RooUnfoldResponse response2D(phi_theta_reco, phi_theta_gen);
	TH1D* X_5D = new TH1D("X_5D", "",nbinx,1, nbinx);
	TH1D* X_5DReco = new TH1D("X_5DReco", "",nbinx,1, nbinx);
	Df_train.Foreach([&response1D, &response2D, &response,phi_theta_data_tr,X_5D,X_5DReco, h_5D_xF, h_5D_pT, h_5D_mass, h_5D_costh, h_5D_phi](bool pass, float reco_mass, float tr_mass, float reco_phi, float tr_phi, float reco_costh, float tr_costh, float reco_xF,  float tr_xF, float reco_pT, float tr_pT, float wt, int trig){ 

			int i_xF= h_5D_xF->GetXaxis()->FindBin(tr_xF);
			int i_xF_= h_5D_xF->GetXaxis()->FindBin(reco_xF);
			int i_pT= h_5D_pT->GetXaxis()->FindBin(tr_pT);
			int i_pT_= h_5D_pT->GetXaxis()->FindBin(reco_pT);
			int i_mass= h_5D_mass->GetXaxis()->FindBin(tr_mass);
			int i_mass_= h_5D_mass->GetXaxis()->FindBin(reco_mass);
			int i_phi= h_5D_phi->GetXaxis()->FindBin(tr_phi);
			int i_phi_= h_5D_phi->GetXaxis()->FindBin(reco_phi);
			int i_costh= h_5D_costh->GetXaxis()->FindBin(tr_costh);
			int i_costh_= h_5D_costh->GetXaxis()->FindBin(reco_costh);

			int n1 = h_5D_xF->GetNbinsX();
			int n2 = h_5D_pT->GetNbinsX();
			int n3 = h_5D_mass->GetNbinsX();
			int n4 = h_5D_phi->GetNbinsX();
			int n5 = h_5D_costh->GetNbinsX();
			int BinIndex1D_true = (i_xF-1) * (n2 * n3 * n4 * n5) + (i_pT-1) * (n3 * n4 * n5) + (i_mass-1) * (n4 * n5) + (i_phi-1) *n5 + (i_costh-1)+1;
			int BinIndex1D_reco = (i_xF_-1) * (n2 * n3 * n4 * n5) + (i_pT_-1) * (n3 * n4 * n5) + (i_mass_-1) * (n4 * n5) + (i_phi_-1) *n5 + (i_costh_-1)+1;

			//for 2D to 1D
			int xT= (( i_phi- 1) * nBinCosth) + i_costh; //true 2D->1D
			int xM =  (( i_phi_- 1) * nBinCosth) + i_costh_; //Measured 2D-->1D

			X_5D->Fill(BinIndex1D_true, wt);
			if(pass==true)X_5DReco->Fill(BinIndex1D_reco, wt);
			if(pass==true) {
				response.Fill(BinIndex1D_reco,BinIndex1D_true,wt);
				response2D.Fill(reco_phi,reco_costh,tr_phi,tr_costh,wt); 
				response1D.Fill(xM, xT,wt); }

			else { 
				response1D.Miss(xT,wt);
				response.Miss(BinIndex1D_true,wt);
				response2D.Miss(tr_phi,tr_costh,wt);
			}
	}, {"passed","mass","true_mass", "phi", "true_phi", "costh", "true_costh", "xF", "true_xF", "pT", "true_pT", "weight", "fpga1"} );
	TFile *acc_file = TFile::Open("testMC2D.root","read");
	TFile *acc_file2d = TFile::Open("testData2D.root","read");


	auto reco_test = (TH2D*)acc_file->Get(Form("costh_test_%i",1));
	RooUnfoldBayes unfoldBayes(&response, reco_test,4);

	TH1D* unfolded_bayes = (TH1D*) unfoldBayes.Hunfold();
	auto* c5 = new TCanvas("c5", "c5", 1800, 400);
	c5->cd();
	X_5D->Draw("Hist");
	X_5D->SetTitle("True X_{1d}; X_{1d}; Count");
	c5->SaveAs("plot/X_5D.png");
	c5->cd();
	X_5DReco->SetTitle("Reconstructed X_{1d}; X_{1d}; Count");
	X_5DReco->Draw("Hist");
	c5->SaveAs("plot/X_5D_reco.png");
	reco_test->Draw("Hist");
	c5->SaveAs("plot/reco_test.png");
	auto* R = response.HresponseNoOverflow();
	c5->cd();
	unfolded_bayes->Draw("hist");
	c5->SaveAs("plot/unfolded_bayes.png");

	TCanvas* c6 = new TCanvas("c6", "c6", 600, 400);
	//R->Draw("colz");

	c6->cd();
	R->SetStats(0);
	R->Draw("colz");
	R->SetTitle("Response Matrix of X_{1d};True X_{1d};Reco X_{1d}");
	c6->SaveAs("plot/testresponse.png");
	c6->cd();

	R->SetTitle("Response Matrix of X_{1d};True X_{1d};Reco X_{1d}");
	R->GetXaxis()->SetRange(1,50);
	R->GetYaxis()->SetRange(1,50);
	R->Draw("colz");
	c6->SaveAs("plot/response.png");


	float nu_high=1.2;
	TH1D* test_lambda_bayes = new TH1D("test_lambda_bayes", "test_lambda_bayes", 20, -0.5, 2.0);
	TH1D* test_mu_bayes = new TH1D("test_mu_bayes", "test_mu_bayes", 20, -0.3, nu_high);
	TH1D* test_nu_bayes = new TH1D("test_nu_bayes", "test_nu_bayes", 20, -0.3, nu_high);


	TH1D* test_lambda_bayes2D = new TH1D("test_lambda_bayes2D", "test_lambda_bayes", 20, -0.5, 2.0);
	TH1D* test_mu_bayes2D = new TH1D("test_mu_bayes2D", "test_mu_bayes", 20, -0.3, nu_high);
	TH1D* test_nu_bayes2D = new TH1D("test_nu_bayes2D", "test_nu_bayes", 20, -0.3, nu_high);

	TH1D* test_lambda_bayes1D = new TH1D("test_lambda_bayes1D", "test_lambda_bayes", 20, -0.5, 2.0);
	TH1D* test_mu_bayes1D = new TH1D("test_mu_bayes1D", "test_mu_bayes", 20, -0.3, nu_high);
	TH1D* test_nu_bayes1D = new TH1D("test_nu_bayes1D", "test_nu_bayes", 20, -0.3, nu_high);

	TF2* fit2D = new TF2("fit2D", "[0] * ( 1 + [1]*y*y + 2*[2]*sqrt(1-y*y)*y*cos(x) + [3]*(1-y*y)*cos(2*x)/2.) ", -M_PI, M_PI,cos_low, cos_high);
	fit2D->SetParNames("A", "#lambda","#mu","#nu");
	fit2D->SetParameters(1,1,0,0);

	TF2* fit2D_2d = new TF2("fit2D_2d", "[0] * ( 1 + [1]*y*y + 2*[2]*sqrt(1-y*y)*y*cos(x) + [3]*(1-y*y)*cos(2*x)/2.) ", -M_PI, M_PI,cos_low, cos_high);
	fit2D_2d->SetParNames("A", "#lambda","#mu","#nu");
	fit2D_2d->SetParameters(1,1,0,0);
	TF2* fit2D_1d = new TF2("fit2D_1d", "[0] * ( 1 + [1]*y*y + 2*[2]*sqrt(1-y*y)*y*cos(x) + [3]*(1-y*y)*cos(2*x)/2.) ", -M_PI, M_PI,cos_low, cos_high);
	fit2D_1d->SetParNames("A", "#lambda","#mu","#nu");
	fit2D_1d->SetParameters(1,1,0,0);
	gStyle->SetOptFit(1);
	bool fitOK=false;
	int nitr =6;
	TH1D* h_reco_test;
	for(int k=1; k<=41; k++){
		h_reco_test = (TH1D*)acc_file->Get(Form("costh_test_%i",k));
		auto reco_test2D = (TH2D*)acc_file2d->Get(Form("costh_test_%i",k));
		TH1D*phi_costh_1d = Make2Dto1D("phi_costh_1d",reco_test2D);
		TH1D* lambda_bayes2d = new TH1D("lambda_bayes2d", "Unfolded #lambda results in different iterations; Iterations; #lambda",nitr,0.5,nitr+0.5);
		TH1D* mu_bayes2d = new TH1D("mu_bayes2d", "Unfolded #mu results in different iterations; Iterations; #mu",nitr,0.5,nitr+0.5);
		TH1D* nu_bayes2d = new TH1D("nu_bayes2d", "Unfolded #nu results in different iterations; Iterations; #nu",nitr,0.5,nitr+0.5);
		TH1D* chi2NDF_bayes2d = new TH1D("chi2NDF_bayes2d","Unfolded #chi^2/NDF results in different iterations; Iterations; #Unfolded #chi^2/NDF",nitr,0.5,nitr+0.5);    

		TH1D* lambda_bayes1d = new TH1D("lambda_bayes1d", "Unfolded #lambda results in different iterations; Iterations; #lambda",nitr,0.5,nitr+0.5);
		TH1D* mu_bayes1d = new TH1D("mu_bayes1d", "Unfolded #mu results in different iterations; Iterations; #mu",nitr,0.5,nitr+0.5);
		TH1D* nu_bayes1d = new TH1D("nu_bayes1d", "Unfolded #nu results in different iterations; Iterations; #nu",nitr,0.5,nitr+0.5);
		TH1D* chi2NDF_bayes1d = new TH1D("chi2NDF_bayes1d","Unfolded #chi^2/NDF results in different iterations; Iterations; #Unfolded #chi^2/NDF",nitr,0.5,nitr+0.5);
		TH1D* lambda_bayes = new TH1D("lambda_bayes", "Unfolded #lambda results in different iterations; Iterations; #lambda",nitr,0.5,nitr+0.5);
		TH1D* mu_bayes = new TH1D("mu_bayes", "Unfolded #mu results in different iterations; Iterations; #mu",nitr,0.5,nitr+0.5);
		TH1D* nu_bayes = new TH1D("nu_bayes", "Unfolded #nu results in different iterations; Iterations; #nu",nitr,0.5,nitr+0.5);
		TH1D* chi2NDF_bayes = new TH1D("chi2NDF_bayes","Unfolded #chi^2/NDF results in different iterations; Iterations; #Unfolded #chi^2/NDF",nitr,0.5,nitr+0.5);  

		float Chisq=1.5;
		int fit_iter =10;

		for(int ii =1; ii<=10; ii++){
			RooUnfoldBayes unfoldBayes(&response, h_reco_test,ii);
			TH1D* unfolded_bayes = (TH1D*) unfoldBayes.Hunfold();
			THnSparseD* h_5D = fill_5D_histogram(unfolded_bayes, nBinxF, xFMin, xFMax, nBinpT, pTMin, pTMax, nBinMass, MassMin, MassMax, nBinPhi, phiMin, phiMax, nBinCosth, costhMin, costhMax);
			//h_5D->GetAxis(3)->SetRange(0,2);
			TH2D* h3proj = h_5D->Projection(4,3); 
			float  tempChisq=0.0;
			int iter=0;
			fitOK=false;
			while(!fitOK)
			{
				iter = iter+1;
				h3proj->Fit("fit2D","R");
				h3proj->Draw("colz");
				tempChisq = fit2D->GetChisquare() / fit2D->GetNDF();
				fitOK = (tempChisq < 2.0 || iter >10 ) ? true : false;
			}
			if(tempChisq < 2.0 && ii==4) {
				test_lambda_bayes->Fill(fit2D->GetParameter(1));
				test_mu_bayes->Fill(fit2D->GetParameter(2));
				test_nu_bayes->Fill(fit2D->GetParameter(3));
			}
			lambda_bayes->SetBinContent(ii,fit2D->GetParameter(1));
			mu_bayes->SetBinContent(ii,fit2D->GetParameter(2));
			nu_bayes->SetBinContent(ii,fit2D->GetParameter(3));
			lambda_bayes->SetBinError(ii,fit2D->GetParError(1));
			mu_bayes->SetBinError(ii,fit2D->GetParError(2));
			nu_bayes->SetBinError(ii,fit2D->GetParError(3));
			chi2NDF_bayes->SetBinContent(ii, fit2D->GetChisquare() / fit2D->GetNDF());
			//Bayes 2D-->1D
			RooUnfoldBayes unfoldBayes1d(&response1D,phi_costh_1d,ii);
			TH1D* unfolded1d_bayes = (TH1D*) unfoldBayes1d.Hunfold();
			TH2D*h_unfoldBayes2d1d = Make1Dto2D("h_unfoldBayes2d1d",unfolded1d_bayes, phi_theta_data_tr);
			tempChisq=0.0;
			iter=0;
			fitOK=false;
			while(!fitOK)
			{
				iter = iter+1;
				h_unfoldBayes2d1d->Fit("fit2D_1d","R");
				tempChisq = fit2D_1d->GetChisquare() / fit2D_1d->GetNDF();
				fitOK = (tempChisq < Chisq || iter >fit_iter ) ? true : false;
			}
			if(tempChisq < 2.0 && ii==4) {
				test_lambda_bayes1D->Fill(fit2D_1d->GetParameter(1));
				test_mu_bayes1D->Fill(fit2D_1d->GetParameter(2));
				test_nu_bayes1D->Fill(fit2D_1d->GetParameter(3));
			}

			lambda_bayes1d->SetBinContent(ii,fit2D_1d->GetParameter(1));
			mu_bayes1d->SetBinContent(ii,fit2D_1d->GetParameter(2));
			nu_bayes1d->SetBinContent(ii,fit2D_1d->GetParameter(3));
			lambda_bayes1d->SetBinError(ii,fit2D_1d->GetParError(1));
			mu_bayes1d->SetBinError(ii,fit2D_1d->GetParError(2));
			nu_bayes1d->SetBinError(ii,fit2D_1d->GetParError(3));
			chi2NDF_bayes1d->SetBinContent(ii, fit2D_1d->GetChisquare() / fit2D_1d->GetNDF());


			//Bayes 2D
			RooUnfoldBayes unfoldBayes2d(&response2D,reco_test2D ,ii);
			TH2D* unfolded2d_bayes = (TH2D*) unfoldBayes2d.Hunfold();
			tempChisq=0.0;
			iter=0;
			fitOK=false;
			while(!fitOK)
			{
				iter = iter+1;
				unfolded2d_bayes->Fit("fit2D_2d","R");
				tempChisq = fit2D_2d->GetChisquare() / fit2D_2d->GetNDF();
				fitOK = (tempChisq < Chisq || iter >fit_iter ) ? true : false;
			}

			if(tempChisq < 2.0 && ii==4) {
				test_lambda_bayes2D->Fill(fit2D_2d->GetParameter(1));
				test_mu_bayes2D->Fill(fit2D_2d->GetParameter(2));
				test_nu_bayes2D->Fill(fit2D_2d->GetParameter(3));
			}
			lambda_bayes2d->SetBinContent(ii,fit2D_2d->GetParameter(1));
			mu_bayes2d->SetBinContent(ii,fit2D_2d->GetParameter(2));
			nu_bayes2d->SetBinContent(ii,fit2D_2d->GetParameter(3));
			lambda_bayes2d->SetBinError(ii,fit2D_2d->GetParError(1));
			mu_bayes2d->SetBinError(ii,fit2D_2d->GetParError(2));
			nu_bayes2d->SetBinError(ii,fit2D_2d->GetParError(3));
			chi2NDF_bayes2d->SetBinContent(ii, fit2D_2d->GetChisquare() / fit2D_2d->GetNDF());

			c6->cd();
			h3proj->SaveAs("plot/h3proj.png");	
			c6->SaveAs(Form("plot/h3proj_%i.png", k));
			delete  unfolded_bayes, h_5D;
			if (gROOT && gROOT->FindObjectAny("unfolded_bayes"))delete gROOT->FindObjectAny("unfolded_bayes");
			delete unfolded2d_bayes;
		}

		TCanvas * c10= new TCanvas("c10","c10", 1200,800);
		c10->Divide(2,2);
		c10->cd(1);

		float lambda_max= 1.2;
		float lambda_min= -0.4;

		lambda_bayes->SetLineColor(kBlue);
		lambda_bayes->Draw("E");
		lambda_bayes->SetLineWidth(3);
		lambda_bayes->SetFillStyle(0);

		lambda_bayes2d->Draw("E2 same");
		lambda_bayes2d->SetLineColor(kAzure+10);
		lambda_bayes2d->SetFillStyle(0);

		lambda_bayes1d->Draw("E same");
		lambda_bayes1d->SetLineColor(kRed);
		lambda_bayes1d->SetMarkerStyle(45);
		c10->cd(2);
		mu_bayes->SetMaximum(0.5);
		mu_bayes->SetMinimum(-0.1);
		mu_bayes->SetLineColor(kBlue);
		mu_bayes->SetLineWidth(3);
		mu_bayes->SetFillStyle(0);

		mu_bayes->Draw("E");
		mu_bayes2d->Draw("E2 same");
		mu_bayes2d->SetLineColor(kAzure+10);
		mu_bayes2d->SetFillStyle(0);
		mu_bayes1d->Draw("E same");
		mu_bayes1d->SetLineColor(kRed);
		mu_bayes1d->SetMarkerStyle(45);
		//setting the legend
		auto leg = new TLegend(0.55,0.92,0.90,0.75);
		leg->AddEntry(mu_bayes,"Bayesian Iterative 5D->1D","lep");
		leg->AddEntry(mu_bayes2d,"Bayesian Iterative 2D","f");
		leg->AddEntry(mu_bayes1d,"Bayesian Iterative 2D->1D","lep");
		leg ->Draw();
		//
		c10->cd(3);
		nu_bayes->SetMinimum(-0.1);
		nu_bayes->SetMaximum(0.4);
		nu_bayes->SetMarkerColor(kBlue);
		nu_bayes->SetLineColor(kBlue);
		nu_bayes->SetLineWidth(3);
		nu_bayes->Draw("E");

		nu_bayes->SetFillStyle(0);
		nu_bayes2d->Draw("E2 same");
		nu_bayes2d->SetLineColor(kAzure+10);
		nu_bayes2d->SetFillStyle(0);

		nu_bayes1d->Draw("E same");
		nu_bayes1d->SetLineColor(kRed);
		nu_bayes1d->SetMarkerStyle(45);

		c10->cd(4);
		chi2NDF_bayes->SetMinimum(0);
		chi2NDF_bayes->Draw("HIST");
		chi2NDF_bayes->SetLineColor(kBlue);

		chi2NDF_bayes2d->SetLineColor(kAzure+10);
		chi2NDF_bayes2d->Draw("HIST same");

		chi2NDF_bayes1d->SetLineColor(kBlack);
		chi2NDF_bayes1d->SetMarkerStyle(41);
		chi2NDF_bayes1d->Draw("HIST same");     

		c10->SaveAs(Form("plot/par_%i.png",k));


	}
	gStyle->SetOptStat(1111);
	TCanvas* c4 = new TCanvas("c4","c4",1800,600);
	c4->Divide(3,1);
	c4->cd(1);
	test_lambda_bayes->Draw("HIST");
	test_lambda_bayes->SetMaximum(35);
	test_lambda_bayes->SetLineColor(kBlue);
	test_lambda_bayes->SetLineWidth(2);
	test_lambda_bayes->SetTitle(";#lambda; Yield");
	c4->cd(2);
	test_mu_bayes->Draw("HIST");
	test_mu_bayes->SetMaximum(35);
	test_mu_bayes->SetLineColor(kBlue);
	test_mu_bayes->SetLineWidth(2);
	test_mu_bayes->SetTitle(";#mu; Yield");
	c4->cd(3);
	test_nu_bayes->Draw("HIST");
	test_nu_bayes->SetMaximum(35);
	test_nu_bayes->SetLineColor(kBlue);
	test_nu_bayes->SetLineWidth(2);
	test_nu_bayes->SetTitle(";#nu; Yield");
	c4->SaveAs("plot/test_bayes.png");
	c4->SaveAs("plot/test_bayes.pdf");

	TCanvas* c15 = new TCanvas("c15","c15",1800,600);
	c15->Divide(3,1);
	c15->cd(1);
	test_lambda_bayes2D->Draw("HIST");
	test_lambda_bayes2D->SetMaximum(35);
	test_lambda_bayes2D->SetLineColor(kBlue);
	test_lambda_bayes2D->SetLineWidth(2);
	test_lambda_bayes2D->SetTitle(";#lambda; Yield");
	c15->cd(2);
	test_mu_bayes2D->Draw("HIST");
	test_mu_bayes2D->SetMaximum(35);
	test_mu_bayes2D->SetLineColor(kBlue);
	test_mu_bayes2D->SetLineWidth(2);
	test_mu_bayes2D->SetTitle(";#mu; Yield");
	c15->cd(3);
	test_nu_bayes2D->Draw("HIST");
	test_nu_bayes2D->SetMaximum(35);
	test_nu_bayes2D->SetLineColor(kBlue);
	test_nu_bayes2D->SetLineWidth(2);
	test_nu_bayes2D->SetTitle(";#nu; Yield");
	c15->SaveAs("plot/test_bayes2D.png");
	c15->SaveAs("plot/test_bayes2D.pdf");

	TCanvas* c16 = new TCanvas("c16","c16",1800,600);
	c16->Divide(3,1);
	c16->cd(1);
	test_lambda_bayes1D->Draw("HIST");
	test_lambda_bayes1D->SetMaximum(35);
	test_lambda_bayes1D->SetLineColor(kBlue);
	test_lambda_bayes1D->SetLineWidth(2);
	test_lambda_bayes1D->SetTitle(";#lambda; Yield");
	c16->cd(2);
	test_mu_bayes1D->Draw("HIST");
	test_mu_bayes1D->SetMaximum(35);
	test_mu_bayes1D->SetLineColor(kBlue);
	test_mu_bayes1D->SetLineWidth(2);
	test_mu_bayes1D->SetTitle(";#mu; Yield");
	c16->cd(3);
	test_nu_bayes1D->Draw("HIST");
	test_nu_bayes1D->SetMaximum(35);
	test_nu_bayes1D->SetLineColor(kBlue);
	test_nu_bayes1D->SetLineWidth(2);
	test_nu_bayes1D->SetTitle(";#nu; Yield");
	c16->SaveAs("plot/test_bayes1D.png");
	c16->SaveAs("plot/test_bayes1D.pdf");
}
