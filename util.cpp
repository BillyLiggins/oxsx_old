#include "util.h"

TH1D* diffHist(TH1D * h1,TH1D * h2){


        double minBin=h1->GetXaxis()->GetXmin();
        double maxBin=h1->GetXaxis()->GetXmax();
        double sliceWidth=h1->GetXaxis()->GetBinWidth(1);
        double numOfBins=h1->GetNbinsX();
	std::cout<<"minBin = "<<minBin<<std::endl;
	std::cout<<"maxBin = "<<maxBin<<std::endl;
	std::cout<<"sliceWidth = "<<sliceWidth<<std::endl;
	
	TH1D* rhist = new TH1D("rhist","",numOfBins,minBin,maxBin);
	for(double i=0;i<numOfBins;i++){
		double h1cont=h1->GetBinContent(i);
		double h2cont=h2->GetBinContent(i);
		double weight= (h1cont-h2cont)/h1cont;
		rhist->SetBinContent(i,weight);
	}

	return rhist;


}

//TH1D doublePlot( const std::vector<TH1D>& histList){




}
