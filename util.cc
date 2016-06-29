#include "util.h"
#include <BinnedPdf.h>        
#include <Rand.h>
#include <string>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <GridSearch.h>
#include <MetropolisHastings.h> 
#include <Minuit.h> 
#include <Histogram.h>
#include <PdfConverter.h>        
#include <CompositePdf.h>        
#include <Convolution.h>        
#include <Gaussian.h>        
#include <DataSetGenerator.h>        
#include <OXSXDataSet.h>        
#include <BoxCut.h>        
#include <vector>        
#include <TH1D.h>        
#include <TH2D.h>        
#include <TFrame.h> 
#include <TGraph.h> 
#include <TGraphErrors.h> 
#include <TPaveText.h> 
#include <TCanvas.h> 
#include <TSystemDirectory.h> 
#include <TList.h> 
#include <THStack.h> 
#include <TLegend.h> 
#include <TLine.h> 
#include <TFile.h> 
#include <TNtuple.h> 
#include <TStyle.h> 
#include <THStack.h> 
#include <TPad.h> 
#include <TPaveStats.h> 
#include <TAttFill.h>
#include <math.h>	
#include <iomanip>


double find2n2bRate(double percentage=0.5,double fidRad=6000)
{
				//This function returns the 2n2b rate given a % loading and fid radius.
				double Vratio= pow(fidRad,3)/pow(6000,3);	
				std::cout << "Volume ratio = "<<Vratio << std::endl;
				double massOfLABPPO=(7.8e5)*Vratio;
				std::cout << "Mass of LABPPO = "<<massOfLABPPO << std::endl;
				double rawTeMass= percentage*massOfLABPPO/100;
				std::cout << "Raw mass of Te = "<<rawTeMass << std::endl;
				double _130TeMass= 0.3408*rawTeMass;
				std::cout << "130Te mass = "<<_130TeMass << std::endl;
				double _130TeHalflife= 7e20;
				std::cout << "130Te halflife = "<<_130TeHalflife << std::endl;
				double decayrate=log(2)/_130TeHalflife;
				std::cout << "decay rate = "<<decayrate << std::endl;
				double molarMassTe=129.906;
				std::cout << "130Te molar mass = "<<molarMassTe << std::endl;
				double molesOfTe=molarMassTe*_130TeMass;
				std::cout << "Number of moles of 130Te = "<<molesOfTe << std::endl;
				double numberOfParticles=molesOfTe*molarMassTe;
				std::cout << "num of particles = "<<numberOfParticles<< std::endl;
				double numOfEvents=numberOfParticles*decayrate*6.022e23;
				std::cout << "Number of events  = "<<numOfEvents  << std::endl;
				std::cout << "check this function again."  << std::endl;

				return numOfEvents;
}

double scaleForTime(double yearRate,double runtime){
				//Should be given the runtime in days and yearly rates (events per year).
				return yearRate*runtime/365;

}

std::vector<double> normRates(std::vector<double>& rates , double normConstant){
				//Should be given the runtime in days and yearly rates (events per year).
				double highestRate=*max_element(rates.begin(),rates.end());
				std::vector<double> scaledRates;
				std::cout << "highest rate = "<<highestRate << std::endl;
				for (int i = 0; i < rates.size(); i++) {
								scaledRates.push_back(rates[i]*normConstant/highestRate);	
				}
				return scaledRates;

}
TH1D* diffHist(TH1D * h1,TH1D * h2){


				double minBin=h1->GetXaxis()->GetXmin();
				double maxBin=h1->GetXaxis()->GetXmax();
				double sliceWidth=h1->GetXaxis()->GetBinWidth(1);
				double numOfBins=h1->GetNbinsX();
				std::cout<<"minBin = "<<minBin<<std::endl;
				std::cout<<"maxBin = "<<maxBin<<std::endl;
				std::cout<<"sliceWidth = "<<sliceWidth<<std::endl;
				std::cout<<"number of bins = "<<numOfBins<<std::endl;

				TH1D* rhist = new TH1D("rhist","",numOfBins,minBin,maxBin);
				for(double i=0;i<numOfBins;i++){
								double h1cont=h1->GetBinContent(i);
								double h2cont=h2->GetBinContent(i);
								double weight;
								if (h1cont!=0 && h1cont-h2cont!=0) {
												weight= (h1cont-h2cont)/h1cont;
								} else {
												weight= 0;
								}
								rhist->SetBinContent(i,weight);
								// std::cout << "weight = "<<weight << std::endl;

				}
				// TCanvas * c1 =new TCanvas();
				// rhist->Draw();
				// c1->Print("temp.png");

				return rhist;


}

void normChecker(std::vector<double>& expectedrate, std::vector<double>& fit_result){
				// std::setw(2);
				// std::setprecision(5);
				std::cout.precision(5);
				std::cout<<"Comparison between expected rates and fit rates:"<<std::endl;
				std::cout<<"------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
				for (int i = 0; i < expectedrate.size(); i++) {
					std::cout<<"Expected rate:\t"	
									<<expectedrate[i]<<"\t|\t fit rate:\t"
									<<fit_result[i]<<"\t|\t Abs difference:\t"
									<<fabs(expectedrate[i]-fit_result[i])<<"\t|\t Frac error:\t"
									<<(expectedrate[i]-fit_result[i])/expectedrate[i]
									<<std::endl;
					std::cout<<"------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
				}
				int numberOfSystmatics= fit_result.size()-expectedrate.size();
				if (numberOfSystmatics){
				}
}
