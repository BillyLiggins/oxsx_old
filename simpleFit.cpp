// A simple fit in energy for signal and a background        
#include <iostream>        
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

using namespace std;

const std::string bgMCfile    = "testData/TeLoadedTe130_2n2b_r2820_s0_p1.ntuple_oxsx.root";
/* const std::string bgMCfile    = "testData/TeLoadedTe130_0n2b_r23_s0_p1.ntuple_oxsx.root"; */
const std::string sigMCfile   = "testData/TeLoadedTe130_0n2b_r23_s0_p1.ntuple_oxsx.root";
const std::string bgTreeName  = "output";
const std::string sigTreeName = "output";

/* const std::string dataFile    = "testData/TeLoadedTe130_0n2b_r23_s0_p1.ntuple_oxsx.root"; */
const std::string dataFile = "combinedDataSet.ntuple_oxsx.root";
const std::string dataTreeName = "output";


vector<std::string> glob( const std::string& path, const std::string& start )
{
  vector<std::string> result;
  TSystemDirectory dir(path.c_str(), path.c_str());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith( start ) ) {
        stringstream fullPath; fullPath << path << "/" << fname.Data();
        result.push_back(fullPath.str());
      }
    }
  }
  return result;
}

void makePDFS(const vector<std::string>& listOfFiles, BinnedPdf & pdf){
	/* for(int j=0;listOfFiles.size();j++){ */
	for(int j=0;j<20;j++){
		cout<<"File number : "<<j<<endl;
		cout<<"File name : "<<listOfFiles[j].c_str()<<endl;

	    ROOTNtuple file(listOfFiles[j], "output");
	    for(size_t i = 0; i <file.GetNEntries(); i++){        
		pdf.Fill(file.GetEntry(i));        
	    }//Loop enteries  
	}//Loop over files

}

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
	const vector<string> signalFileList= glob("/data/snoplus/liggins/year1/oxsx/testData","TeLoadedTe130_0n2b_r"); 
	const vector<std::string> backgroundFileList= glob("/data/snoplus/liggins/year1/oxsx/testData","TeLoadedTe130_2n2b_r"); 
	 
	cout<<"File name after glob : "<<signalFileList[0]<<endl;
	makePDFS(signalFileList,signalPdf);
	makePDFS(backgroundFileList,bgPdf);

	/* ROOTNtuple bgMC(bgMCfile, bgTreeName); */
	/* ROOTNtuple signalMC(sigMCfile, sigTreeName); */
	 
	/* for(size_t i = 0; i < bgMC.GetNEntries(); i++){ */        
	/*     bgPdf.Fill(bgMC.GetEntry(i)); */        
	/* } */        
	  
	/* for(size_t i = 0; i < signalMC.GetNEntries(); i++){ */        
	/*     signalPdf.Fill(signalMC.GetEntry(i)); */        
	/* } */        
	     
	 
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

	std::cout << "Start Search"<<std::endl; 
	FitResult result = gSearch.Optimise(&lhFunction);
	std::cout << "End Search"<<std::endl; 
	 
	std::vector<double> fit = result.GetBestFit();        
	result.Print();
	result.SaveAs("simpleFit_result.txt");
	return 0;        
}
