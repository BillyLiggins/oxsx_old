// A simple fit in energy for signal and a background        
#include <BinnedPdf.h>        
#include <string.h>        
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
#include<BoxCut.h>        

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
#include <TPad.h> 
#include <TPaveStats.h> 
#include <TAttFill.h>

const std::string bgMCfile    = "testData/complete2b2n.ntuple_oxsx.root";
const std::string sigMCfile   = "testData/complete0b2n.ntuple_oxsx.root"; 
const std::string bgTreeName  = "output";
const std::string sigTreeName = "output";

/* const std::string dataFile = "oneFileToRuleThemAll.ntuple.root"; */
/* const std::string dataFile = "secondCompleteData.ntuple.root"; */
const std::string dataFile = "testData/tester.ntuple.root";
/* const std::string dataFile = "full.ntuple.root"; */
/* const std::string dataFile = "secondCompleteData_2.ntuple.root"; */
/* const std::string dataFile = "combinedDataSet.ntuple_oxsx.root"; */
/* const std::string dataFile = "complete2b2n.ntuple_oxsx.root"; */
const std::string dataTreeName = "output";

TH1D* diffHist(TH1D * h1,TH1D * h2){


        double minBin=h1->GetXaxis()->GetXmin();
        double maxBin=h1->GetXaxis()->GetXmax();
        double sliceWidth=h1->GetXaxis()->GetBinWidth(1);
        double numOfBins=h1->GetNbinsX();
	// std::cout<<"minBin = "<<minBin<<std::endl;
	// std::cout<<"maxBin = "<<maxBin<<std::endl;
	// std::cout<<"sliceWidth = "<<sliceWidth<<std::endl;
	
	TH1D* rhist = new TH1D("rhist","",numOfBins,minBin,maxBin);
	for(double i=0;i<numOfBins;i++){
		double h1cont=h1->GetBinContent(i);
		double h2cont=h2->GetBinContent(i);
		double weight= (h1cont-h2cont)/h1cont;
		rhist->SetBinContent(i,weight);
	}

	return rhist;


}


int main(){        
	////////////////////        1. Set Up PDFs //        //////////////////        

	// Set up binning        
	AxisCollection axes;
	axes.AddAxis(PdfAxis("energy", 0., 3, 60,"Energy"));
	// axes.AddAxis(PdfAxis("energy", 2., 3, 20,"Energy"));

	// Only interested in first bit of data ntuple        
	DataRepresentation dataRep(0);        

	// Set up pdf with these bins in this observable        
	BinnedPdf bgPdf(axes);
	bgPdf.SetDataRep(dataRep);

	BinnedPdf signalPdf(axes);
	signalPdf.SetDataRep(dataRep);

	BinnedPdf comPdf(axes);
	comPdf.SetDataRep(dataRep);        

	std::cout << "Initialised Pdfs" << std::endl;        

	/////////////////////////////////////        2. Fill with data and
	//normalise //        ///////////////////////////////////        

	ROOTNtuple bgMC(bgMCfile, bgTreeName); ROOTNtuple signalMC(sigMCfile,sigTreeName);

	BoxCut boxCut(0,0,3); for(size_t i = 0; i < bgMC.GetNEntries(); i++){
		if(boxCut.PassesCut(bgMC.GetEntry(i))){
			bgPdf.Fill(bgMC.GetEntry(i));        
		}
	}        

	for(size_t i = 0; i < signalMC.GetNEntries(); i++){
		if(boxCut.PassesCut(signalMC.GetEntry(i))){
			signalPdf.Fill(signalMC.GetEntry(i));        
		}
	}        


	bgPdf.Normalise();
	signalPdf.Normalise();        

	TCanvas * sigCan = new TCanvas();
	TH1D  sigPlot =	PdfConverter::ToTH1D(signalPdf,false);
	sigPlot.Draw();
	sigCan->Print("sigPlot.png");

	TCanvas * bgCan = new TCanvas();
	TH1D  bgPlot =PdfConverter::ToTH1D(bgPdf,false);
	bgPlot.Draw();
	bgCan->Print("bgPlot.png");


	std::cout << "Filled pdfs " << std::endl;        

	//////////////////////////// 3. Set Up LH function //////////////////////////////        

	ROOTNtuple dataNt(dataFile, dataTreeName);

	for(size_t i = 0; i < dataNt.GetNEntries(); i++){        
		//Cut between 0 and 3.
		if(boxCut.PassesCut(dataNt.GetEntry(i))){
			comPdf.Fill(dataNt.GetEntry(i));        
		}
	}        

	std::cout<< "n = " << comPdf.Integral()<< std::endl;
	/* comPdf.Normalise(); */

	/* return 0; */
	TCanvas * comCan = new TCanvas();
	TH1D  comPlot =PdfConverter::ToTH1D(comPdf,false); 
	comPlot.Draw();
	comCan->Print("comPlot.png");

	TCanvas * normCan = new TCanvas();
	TH1D * comPlotcopy = (TH1D*) comPlot.Clone();
	comPlotcopy->SetLineColor(kGreen);
	comPlotcopy->SetFillColor(kGreen);
	TAxis * xaxis= comPlotcopy->GetXaxis();
	/* TAxis * yaxis= comPlotcopy->GetYaxis(); */
	double scale= comPlotcopy->Integral(xaxis->FindBin(0.),xaxis->FindBin(10.));
	comPlotcopy->Scale(1/scale);
	comPlotcopy->SetMaximum(0.2);
	/* comPlotcopy->DrawNormalized(); */
	comPlotcopy->Draw();
	sigPlot.SetLineColor(kRed);
	sigPlot.SetLineWidth(2);
	sigPlot.Draw("same");
	bgPlot.SetLineColor(kBlue);
	bgPlot.SetLineWidth(2);
	bgPlot.Draw("same");

	TLegend* normleg =new TLegend(0.1,0.7,0.48,0.9);
	normleg->AddEntry(comPlotcopy,"Normalised Data","lf");
	normleg->AddEntry(&sigPlot,"Signal Pdf","lf");
	normleg->AddEntry(&bgPlot,"Background Pdf","lf");
	normleg->Draw();
	normCan->Print("normPlot.png");

	//----
	//

	Convolution conv;
	Gaussian gaus(0,1);
	conv.SetFunction(&gaus);
	conv.SetAxes(axes);
	DataRepresentation DataRep(0);

	conv.SetDataRep(DataRep);
	conv.SetPdfDataRep(DataRep);

	/////////////////// Setting up the lhFunctions ////////////////////////
	BinnedNLLH lhFunction_gSearch;
	/* lhFunction_gSearch.SetDataSet(&dataNt); // initialise withe the data set */
	lhFunction_gSearch.SetBufferAsOverflow(false);        
	lhFunction_gSearch.SetBuffer(0,10,10);
	lhFunction_gSearch.SetDataPdf(comPdf);
	lhFunction_gSearch.AddPdf(bgPdf);
	lhFunction_gSearch.AddPdf(signalPdf);        
	lhFunction_gSearch.AddSystematic(&conv);        

	BinnedNLLH lhFunction_minuit;
	lhFunction_minuit.SetDataSet(&dataNt); // initialise withe the data set
	lhFunction_minuit.AddPdf(bgPdf);
	lhFunction_minuit.AddPdf(signalPdf);        

	BinnedNLLH lhFunction_metHast;
	lhFunction_metHast.SetDataSet(&dataNt); // initialise withe the data set
	lhFunction_metHast.AddPdf(bgPdf);
	lhFunction_metHast.AddPdf(signalPdf);        
	std::cout << "Built LH functions " << std::endl;        


	/////////////////// Set up the optimisations ////////////////////////////
	GridSearch gSearch;        

	std::vector<double> minima;
	/* minima.push_back(90000); */
	/* minima.push_back(90000); */
	minima.push_back(70000);
	minima.push_back(70000);
	minima.push_back(-0.002);
	minima.push_back(0.);
	std::vector<double> maxima;
	maxima.push_back(140000);
	maxima.push_back(140000);
	maxima.push_back(0.0002);
	maxima.push_back(0.5);
	std::vector<double> stepsizes;
	stepsizes.push_back(1);
	stepsizes.push_back(1);
	stepsizes.push_back(0.001);
	stepsizes.push_back(0.001);

	gSearch.SetMaxima(maxima);
	gSearch.SetMinima(minima);
	gSearch.SetStepSizes(stepsizes);        

	Minuit minuit;        
	std::vector<double> InitialValues;
	InitialValues.push_back(100000);
	InitialValues.push_back(100000);
	InitialValues.push_back(0);
	InitialValues.push_back(0.2);

	std::vector<double> InitialErrors;
	InitialErrors.push_back(1);
	InitialErrors.push_back(1);
	InitialErrors.push_back(0.001);
	InitialErrors.push_back(0.001);
	minuit.SetMaxima(maxima);
	minuit.SetMinima(minima);        
	minuit.SetInitialValues(InitialValues);        
	minuit.SetInitialErrors(InitialErrors);        

	MetropolisHastings metHast;        
	/* metHast.SetMaxIter(10000000); */
	metHast.SetMaxIter(10000000);
	metHast.SetMaxima(maxima);
	metHast.SetMinima(minima);        
	metHast.SetFlipSign(true);
	metHast.SetTestStatLogged(true);

	std::vector<double> sigmas(4,10);
	metHast.SetSigmas(sigmas);


	////////////     Now perform the fits ////////////        
	/* FitResult result_gSearch = gSearch.Optimise(&lhFunction_gSearch); */
	FitResult result_minuit = minuit.Optimise(&lhFunction_gSearch);
	FitResult result_metHast = metHast.Optimise(&lhFunction_gSearch);

	/* std::vector<double> fit_gSearch = result_gSearch.GetBestFit(); */        
	std::vector<double> fit_gSearch(2,1); 
	/* result_gSearch.Print(); */
	/* result_gSearch.SaveAs("simpleFit_result_gSearch.txt"); */

	std::vector<double> fit_minuit = result_minuit.GetBestFit();        
	std::cout<<"Minuit Results, pdf 0 = bkg, 1 = sig"<<std::endl;
	result_minuit.Print();
	result_minuit.SaveAs("simpleFit_result_minuit.txt");

	std::vector<double> fit_metHast =result_metHast.GetBestFit();        
	std::cout<<"MetHast Results, pdf 0 = bkg, 1 = sig"<<std::endl;
	/* std::vector<double> fit_metHast(2,1); */ 
	Histogram hist = result_metHast.GetStatSpace();
	/* PdfConverter::ToTH2D(hist).SaveAs("lh_2d.root"); */
	PdfConverter::ToTH1D(hist.Marginalise(0)).SaveAs("pdf_norm0.root");
	PdfConverter::ToTH1D(hist.Marginalise(1)).SaveAs("pdf_norm1.root");
	PdfConverter::ToTH1D(hist.Marginalise(2)).SaveAs("pdf_norm2.root");
	PdfConverter::ToTH1D(hist.Marginalise(3)).SaveAs("pdf_norm3.root");
	result_metHast.Print();
	result_metHast.SaveAs("simpleFit_result_metHast.txt");



	///////////////Setting up the histograms /////////////////////

	/* gStyle->SetOptStat(kFALSE); */ 
	TH1D * scaledBg_gSearch= (TH1D*) bgPlot.Clone();
	TH1D * scaledSig_gSearch= (TH1D*) sigPlot.Clone();
	TH1D * scaledBg_minuit= (TH1D*) bgPlot.Clone();
	TH1D * scaledSig_minuit= (TH1D*) sigPlot.Clone();
	TH1D * scaledBg_metHast= (TH1D*) bgPlot.Clone();
	TH1D * scaledSig_metHast= (TH1D*) sigPlot.Clone();

	scaledBg_gSearch->Scale(fit_gSearch[0]);
	scaledSig_gSearch->Scale(fit_gSearch[1]);

	scaledBg_minuit->Scale(fit_minuit[0]);
	scaledSig_minuit->Scale(fit_minuit[1]);

	scaledBg_metHast->Scale(fit_metHast[0]);
	scaledSig_metHast->Scale(fit_metHast[1]);

	TH1D * complete_gSearch= (TH1D*) scaledBg_gSearch->Clone();
	TH1D * complete_minuit= (TH1D*) scaledBg_minuit->Clone();
	TH1D * complete_metHast=(TH1D*) scaledBg_metHast->Clone();

	complete_gSearch->Add(scaledBg_gSearch,scaledSig_gSearch,1,1);
	complete_minuit->Add(scaledBg_minuit,scaledSig_minuit,1,1);
	complete_metHast->Add(scaledBg_metHast,scaledSig_metHast,1,1);

	TCanvas * fitCan= new TCanvas();

	comPlot.SetFillColorAlpha(kGreen,0.5);
	comPlot.Draw();

	complete_gSearch->SetLineColor(kBlue);
	complete_gSearch->SetLineWidth(2);
	complete_minuit->SetLineColor(kRed);
	complete_minuit->SetLineWidth(2);
	complete_metHast->SetLineColor(kBlack);
	complete_metHast->SetLineWidth(2);

	/* complete_gSearch->SetFillColorAlpha(kBlue,0.5); */
	/* complete_minuit->SetFillColorAlpha(kRed,0.5); */
	/* complete_metHast->SetFillColorAlpha(kBlack,0.5); */

	complete_gSearch->Draw("same");
	complete_minuit->Draw("same");
	complete_metHast->Draw("same");

	TLegend* leg =new TLegend(0.1,0.7,0.48,0.9);
	leg->AddEntry(complete_gSearch,"Grid Search","lf");
	leg->AddEntry(complete_minuit,"Minuit","lf");
	leg->AddEntry(complete_metHast,"Metropolis Hastings","lf");
	leg->Draw();

	fitCan->Print("All_Fits.png");


	//_+_+_+_+_+_+_+_+_+_Double Panel Plots_+_+_+_+_+_+_+_+_+_

	// ------- Create a canvas:
	TCanvas* diff = new TCanvas("diff","",800,800); diff->cd();
	// -------------- Top panel
	gStyle->SetOptStat(kFALSE); 
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	pad1->Draw();
	pad1->cd();
	pad1->SetGrid(kTRUE);
	pad1->SetBottomMargin(0.00);
	gPad->RedrawAxis(); 

	comPlot.Draw();
	/* complete_gSearch->Draw("same"); */
	complete_minuit->Draw("same");
	complete_metHast->Draw("same");
	leg->Draw();
	/* complete_gSearch->GetXaxis()->SetTitle("Energy (MeV)"); */
	/* complete_gSearch->GetYaxis()->SetTitle("Frac error"); */
	diff->cd();
	// -------------- Bottom panel
	TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3);
	pad2->SetTopMargin(0.00);
	pad2->Draw();
	pad2->cd();
	pad2->SetBottomMargin(0.3);
	pad2->SetGrid(kTRUE);
	gStyle->SetOptStat(kFALSE); 

	/* diffHist->Draw(); */

	TH1D * diff_gSearch= diffHist(&comPlot,complete_gSearch);
	TH1D * diff_minuit= diffHist(&comPlot,complete_minuit);
	TH1D * diff_metHast= diffHist(&comPlot,complete_metHast);
	diff_gSearch->SetLineColor(kBlue);
	diff_gSearch->SetLineWidth(2);
	diff_minuit->SetLineColor(kRed);
	diff_minuit->SetLineWidth(2);
	diff_metHast->SetLineColor(kBlack);
	diff_metHast->SetLineWidth(2);

	diff_gSearch->SetFillColorAlpha(kBlue,0.1);
	diff_minuit->SetFillColorAlpha(kRed,0.1);
	diff_metHast->SetFillColorAlpha(kBlack,0.1);

	diff_minuit->GetXaxis()->SetTitle("Energy (MeV)");
	diff_minuit->GetYaxis()->SetTitle("Frac error");
	diff_minuit->GetXaxis()->SetLabelSize(0.1);
	diff_minuit->GetXaxis()->SetTitleSize(0.1);
	diff_minuit->GetYaxis()->SetLabelSize(0.05);
	diff_minuit->GetYaxis()->SetTitleSize(0.1);

	/* diff_gSearch->GetXaxis()->SetTitle("Energy (MeV)"); */
	/* diff_gSearch->GetYaxis()->SetTitle("Frac error"); */
	/* diff_gSearch->GetXaxis()->SetLabelSize(0.1); */
	/* diff_gSearch->GetXaxis()->SetTitleSize(0.1); */
	/* diff_gSearch->GetYaxis()->SetLabelSize(0.1); */
	/* diff_gSearch->GetYaxis()->SetTitleSize(0.1); */

	/* complete_minuit->Draw(); */
	/* diff_gSearch->Draw(); */
	/* diff_minuit->Draw("same"); */
	diff_minuit->Draw();
	/* diff_metHast->Draw("same"); */
	/* gStyle->SetOptStat(0); */ 

	diff->Print("doublePlanel.png");
	std::cout<<" Number diff = "<< comPdf.Integral() - fit_minuit[0]-fit_minuit[1] <<std::endl;

	return 0;
}
