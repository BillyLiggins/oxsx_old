// A simple fit in energy for signal and a background        
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <GridSearch.h>        
#include <Histogram.h>        
#include <PdfConverter.h>        

#include <TH1D.h>        
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
#include <TPad.h>

const std::string bgMCfile    = "complete2b2n.ntuple_oxsx.root";
const std::string sigMCfile   = "complete0b2n.ntuple_oxsx.root";
const std::string bgTreeName  = "output";
const std::string sigTreeName = "output";

const std::string dataFile = "combinedDataSet.ntuple_oxsx.root";
const std::string dataTreeName = "output";

int main(){        
    ////////////////////        
    // 1. Set Up PDFs //        
    ////////////////////        
         
    // Set up binning        
    AxisCollection axes;        
    axes.AddAxis(PdfAxis("energy", 2, 3, 10, "Energy"));
         
    // Only interested in first bit of data ntuple        
    DataRepresentation dataRep(0);        
         
    // Set up pdf with these bins in this observable        
    BinnedPdf bgPdf(axes);      bgPdf.SetDataRep(dataRep);        
    BinnedPdf signalPdf(axes);  signalPdf.SetDataRep(dataRep);        
    BinnedPdf comPdf(axes);  comPdf.SetDataRep(dataRep);        
         
    std::cout << "Initialised Pdfs" << std::endl;        
         
    /////////////////////////////////////        
    // 2. Fill with data and normalise //        
    /////////////////////////////////////        
         
    ROOTNtuple bgMC(bgMCfile, bgTreeName);
    ROOTNtuple signalMC(sigMCfile, sigTreeName);
         
    for(size_t i = 0; i < bgMC.GetNEntries(); i++){        
        bgPdf.Fill(bgMC.GetEntry(i));        
    }        
          
    for(size_t i = 0; i < signalMC.GetNEntries(); i++){        
        signalPdf.Fill(signalMC.GetEntry(i));        
    }        
             
         
    bgPdf.Normalise();        
    signalPdf.Normalise();        
         
	TCanvas * sigCan = new TCanvas();
	TH1D  sigPlot = PdfConverter::ToTH1D(signalPdf,false);
	sigPlot.Draw();
	sigCan->Print("sigPlot.png");
	 
	TCanvas * bgCan = new TCanvas();
	TH1D  bgPlot = PdfConverter::ToTH1D(bgPdf,false);
	bgPlot.Draw();
	bgCan->Print("bgPlot.png");


	std::cout << "Filled pdfs " << std::endl;        
	 
	////////////////////////////        
	// 3. Set Up LH function  //        
	////////////////////////////        
	ROOTNtuple dataNt(dataFile, dataTreeName);


	for(size_t i = 0; i < dataNt.GetNEntries(); i++){        
	comPdf.Fill(dataNt.GetEntry(i));        
	}        
	comPdf.Normalise();

	TCanvas * comCan = new TCanvas();
	TH1D  comPlot = PdfConverter::ToTH1D(comPdf,false);
	comPlot.Draw();
	comCan->Print("comPlot.png");

    BinnedNLLH lhFunction;        
    lhFunction.SetDataSet(&dataNt); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);        
    lhFunction.AddPdf(signalPdf);        
        
    std::cout << "Built LH function " << std::endl;        
         
    // Set up the optimisation        
    GridSearch gSearch;        
             
    std::vector<double> minima;
    minima.push_back(0);
    minima.push_back(0);
    std::vector<double> maxima;
    maxima.push_back(1000);
    maxima.push_back(1000);
    std::vector<double> stepsizes(2, 1);
         
    gSearch.SetMaxima(maxima);        
    gSearch.SetMinima(minima);        
    gSearch.SetStepSizes(stepsizes);        
             
    ////////////        
    // 4. Fit //        
    ////////////        
    FitResult result = gSearch.Optimise(&lhFunction);
         
    std::vector<double> fit = result.GetBestFit();        
    result.Print();
    result.SaveAs("simpleFit_result.txt");
    return 0;        
}
