#ifndef HISTOGRAM_UTILS_H
#define HISTOGRAM_UTILS_H

#include <TH1D.h>
#include <THnSparse.h>

THnSparseD* fill_5D_histogram(TH1D* unfolded_histogram, int nBinxF, double xFMin, double xFMax, int nBinpT, double pTMin, double pTMax, int nBinMass, double MassMin, double MassMax, int nBinPhi, double phiMin, double phiMax, int nBinCosth, double costhMin, double costhMax) {

    // Create the THnSparseD histogram
    const int ndim = 5;
    const int nbins[ndim] = {nBinxF, nBinpT, nBinMass, nBinPhi, nBinCosth};
    const double xmin[ndim] = {xFMin, pTMin, MassMin, phiMin, costhMin};
    const double xmax[ndim] = {xFMax, pTMax, MassMax, phiMax, costhMax};
    THnSparseD* h_5D = new THnSparseD("h_5D", "5D histogram", ndim, nbins, xmin, xmax);
    h_5D->Sumw2();

    // Loop over all bins in the unfolded histogram and fill the corresponding bin in the 5D histogram
    int n1 = nBinxF;
    int n2 = nBinpT;
    int n3 = nBinMass;
    int n4 = nBinPhi;
    int n5 = nBinCosth;

    for (int ixF = 1; ixF <= n1; ixF++) {
        for (int ipT = 1; ipT <= n2; ipT++) {
            for (int imass = 1; imass <= n3; imass++) {
                for (int iphi = 1; iphi <= n4; iphi++) {
                    for (int icosth = 1; icosth <= n5; icosth++) {

                        int x_bin = (ixF-1) * (n2 * n3 * n4 * n5) + (ipT-1) * (n3 * n4 * n5) + (imass-1) * (n4 * n5) + (iphi-1) * n5 + (icosth-1) + 1;
                        int coord[5] = {ixF, ipT, imass, iphi, icosth};
                        int new_content = unfolded_histogram->GetBinContent(x_bin);
                        int new_content_err = unfolded_histogram->GetBinError(x_bin);
                        int bin_idx = h_5D->GetBin(coord);
                        h_5D->SetBinContent(bin_idx, new_content);
                        h_5D->SetBinError(bin_idx, new_content_err);
                    }
                }
            }
        }
    }
    return h_5D;
}

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

#endif // HISTOGRAM_UTILS_H

