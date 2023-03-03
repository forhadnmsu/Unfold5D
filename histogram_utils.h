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

#endif // HISTOGRAM_UTILS_H

