// A simple fit in energy for signal and a background        
#include <iostream>        
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <GridSearch.h>        
#include <Histogram.h>        
#include <PdfConverter.h>        
#include <TH1D.h>        
#include <TCanvas.h>        

const std::string bgMCfile    = "testData/TeLoadedTe130_2n2b_r2820_s0_p1.ntuple_oxsx.root";
const std::string bgTreeName  = "output";

/* const std::string dataFile    = "testData/TeLoadedTe130_0n2b_r23_s0_p1.ntuple_oxsx.root"; */
const std::string dataFile = "combinedDataSet.ntuple_oxsx.root";
/* const std::string dataFile = "testData/TeLoadedTe130_2n2b_r2815_s0_p1.ntuple_oxsx.root"; */
const std::string dataTreeName = "output";

int main(){        
    ////////////////////        
    // 1. Set Up PDFs //        
    ////////////////////        
         
    // Set up binning        
    AxisCollection axes;        
    axes.AddAxis(PdfAxis("energy", 0, 3.0, 30, "Energy"));
         
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
         
    for(size_t i = 0; i < bgMC.GetNEntries(); i++){        
        bgPdf.Fill(bgMC.GetEntry(i));        
    }        
         
    bgPdf.Normalise();        

    /* TCanvas * bgCan = new TCanvas(); */
    /* TH1D  bgPlot = PdfConverter::ToTH1D(bgPdf,false); */
    /* bgPlot.Draw(); */
    /* bgCan->Print("bgPlot.png"); */


    std::cout << "Filled pdfs " << std::endl;        
         
    ////////////////////////////        
    // 3. Set Up LH function  //        
    ////////////////////////////        
    ROOTNtuple dataNt(dataFile, dataTreeName);


    for(size_t i = 0; i < dataNt.GetNEntries(); i++){        
        comPdf.Fill(dataNt.GetEntry(i));        
    }        
    /* comPdf.Normalise(); */

    /* TCanvas * comCan = new TCanvas(); */
    /* TH1D  comPlot = PdfConverter::ToTH1D(comPdf,false); */
    /* comPlot.Draw(); */
    /* comCan->Print("comPlot.png"); */

    BinnedNLLH lhFunction;        
    lhFunction.SetDataSet(&dataNt); // initialise withe the data set
    lhFunction.AddPdf(bgPdf);        
        
    std::cout << "Built LH function " << std::endl;        
         
    // Set up the optimisation        
    GridSearch gSearch;        
             
    std::vector<double> minima;
    minima.push_back(0);
    /* minima.push_back(0); */
    std::vector<double> maxima;
    maxima.push_back(1000);
    /* maxima.push_back(1000); */
    std::vector<double> stepsizes(1, 1);
         
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
    result.SaveAs("signalFit_result.txt");
    return 0;        
}
