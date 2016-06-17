// // A simple fit in energy for 2n2b, Po210 and Bi210 backgrounds.        
// This is very messy you should make it so that the everything is in vectors.

#include "util.h"
#include <BinnedPdf.h>        
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

//const std::string Bi210File    = "completeBi210.ntuple_oxsx.root";
//const std::string Po210File    = "completePo210.ntuple_oxsx.root";
//const std::string DoubleNueFile   = "complete2b2n.ntuple_oxsx.root"; 

//const std::string DoubleNueTreeName  = "output";
//const std::string Bi210TreeName  = "output";
//const std::string Po210TreeName   = "output";

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
//

int main(){        
				bool QmetHast=true;
				std::vector<std::string> inputFiles;
				inputFiles.push_back("completeBi210.ntuple_oxsx.root");
				inputFiles.push_back("completePo210.ntuple_oxsx.root");
				inputFiles.push_back("complete2b2n.ntuple_oxsx.root");
				inputFiles.push_back("testData/TeLoaded/C14/complete_C14.ntuple_oxsx.root");

				std::vector<std::string> names;
				names.push_back("data");
				names.push_back("Bi210");
				names.push_back("Po210");
				names.push_back("2n2b");
				names.push_back("C14");

				std::vector<std::string> inputTrees;
				inputTrees.push_back("output");
				inputTrees.push_back("output");
				inputTrees.push_back("output");
				inputTrees.push_back("output");

				////////////////////        1. Set Up PDFs //      //////////////////        

				// Set up binning        
				AxisCollection axes;
				axes.AddAxis(PdfAxis("energy", 0., 3, 60,"Energy"));

				// Only interested in first bit of data ntuple        
				DataRepresentation dataRep(0);        


				// Set up pdf with these bins in this observable        


				/******************************************************************************************
				 *************************** 2. Fill with data and normalise ******************************
				 ******************************************************************************************/        

				std::vector<double> rates;
				rates.push_back(1000); //Bi210
				rates.push_back(1000); //Po210
				rates.push_back(1000); //2n2b
				rates.push_back(1000); //2n2b
				// rates.push_back(427368); //Bi210
				// rates.push_back(1.74e7); //Po210
				// rates.push_back(100000); //2n2b

				std::vector<ROOTNtuple*> ntupleList;
				std::vector<BinnedPdf> binnedPDFList;
				std::vector<TH1D> histList;
				DataSetGenerator dataGen;

				BoxCut boxCut(0,0,3);
				for(int i =0;i<inputFiles.size();i++){
								std::cout<<"i = "<<i<<std::endl;
								BinnedPdf Pdf(axes);
								std::cout<<"After PDF"<<std::endl;
								binnedPDFList.push_back(Pdf);
								std::cout << binnedPDFList.size()<<" i = "<<i << std::endl;
								binnedPDFList.at(i).SetDataRep(dataRep);


								// ROOTNtuple file_temp(inputFiles.at(1),inputTrees.at(1));
								ntupleList.push_back(new ROOTNtuple(inputFiles.at(i),inputTrees.at(i)));

								dataGen.AddDataSet(ntupleList.at(i),rates.at(i));

								for(size_t j = 0; j < ntupleList.at(i)->GetNEntries(); j++){
												if(boxCut.PassesCut(ntupleList.at(i)->GetEntry(j))){
																binnedPDFList.at(i).Fill(ntupleList.at(i)->GetEntry(j));        
												}
								}        

								binnedPDFList.at(i).Normalise();

								std::cout << "Initialised Pdfs" << std::endl;        

								TCanvas * defCan = new TCanvas();
								histList.push_back(PdfConverter::ToTH1D(binnedPDFList[i],false));
								histList[i].Draw();
								defCan->Print(Form("Pdf_%d.png",i));
				}

				BinnedPdf DataPdf(axes);
				DataPdf.SetDataRep(dataRep);        
				OXSXDataSet fakeData= dataGen.ExpectedRatesDataSet();
				for(size_t i = 0; i < fakeData.GetNEntries(); i++){        
								//Cut between 0 and 3.
								if(boxCut.PassesCut(fakeData.GetEntry(i))){
												DataPdf.Fill(fakeData.GetEntry(i));        
								}
				}        


				TCanvas * comCan = new TCanvas();
				histList.insert(histList.begin(), PdfConverter::ToTH1D(DataPdf,false));
				histList[0].Draw();
				comCan->Print("DataPlot.png");


				TCanvas * pdfCanvas= new TCanvas();
				pdfCanvas->Divide(2,2);
				for(int i=0;i<histList.size();i++){
								pdfCanvas->cd(i+1);
								histList[i].Draw();
				}
				pdfCanvas->Print("allPdfs.png");


				TCanvas * normCan = new TCanvas();
				TH1D * comPlotcopy = (TH1D*) histList[0].Clone();
				comPlotcopy->SetLineColor(kGreen);
				comPlotcopy->SetFillColor(kGreen);
				TAxis * xaxis= comPlotcopy->GetXaxis();
				TAxis * yaxis= comPlotcopy->GetYaxis();
				double scale= comPlotcopy->Integral(xaxis->FindBin(0.),xaxis->FindBin(10.));
				comPlotcopy->Scale(1/scale);
				comPlotcopy->SetMaximum(1);
				// comPlotcopy->DrawNormalized();
				comPlotcopy->Draw();
				TLegend* normleg =new TLegend(0.1,0.7,0.48,0.9);
				normleg->AddEntry(comPlotcopy,"Normalised Data","lf");
				//Start at i=1 because histList[0] is data.

				for (int i = 1; i < histList.size(); i++) {
								std::cout << histList.size() << std::endl;
								histList[i].SetLineColor(i);
								histList[i].SetLineWidth(2);
								histList[i].Draw("same");
								normleg->AddEntry(&histList[i],(names[i]+" pdf").c_str(),"lf");
				}
				normleg->Draw();
				normCan->Print("normPlot.png");



				// Convolution conv;
				// Gaussian gaus(0,1); 
				// conv.SetFunction(&gaus);
				// conv.SetAxes(axes);
				// DataRepresentation DataRep(0);
				// conv.SetDataRep(DataRep);
				// conv.SetPdfDataRep(DataRep);

				///////////////////// Setting up the lhFunctions ////////////////////////

				BinnedNLLH lhFunction_gSearch;
				lhFunction_gSearch.SetDataPdf(DataPdf); // initialise withe the data set
				for(int i=0;i<binnedPDFList.size();i++){
								lhFunction_gSearch.AddPdf(binnedPDFList[i]);

				}
				// lhFunction_gSearch.AddSystematic(&conv);

				// lhFunction_gSearch.SetBufferAsOverflow(false);
				// lhFunction_gSearch.SetBuffer(0,10,10);
				//lhFunction_gSearch.SetDataPdf(DataPdf);
				//lhFunction_gSearch.AddPdf(Bi210Pdf);
				//lhFunction_gSearch.AddPdf(Po210Pdf);
				//lhFunction_gSearch.AddPdf(DoubleNuePdf);
				// lhFunction_gSearch.AddPdf(signalPdf);
				// lhFunction_gSearch.AddSystematic(&conv);
				//
				std::cout << "Built LH functions " << std::endl;


				/////////////////// Set up the optimisations ////////////////////////////
				std::vector<double> minima;
				std::vector<double> maxima;
				std::vector<double> stepsizes;
				std::vector<double> InitialValues;
				std::vector<double> InitialErrors;
				for(int i=0;i<binnedPDFList.size();i++){
								minima.push_back(0);
								maxima.push_back(140000);
								stepsizes.push_back(1);
								InitialValues.push_back(100000);
								InitialErrors.push_back(1);
				}

				Minuit minuit;

				minuit.SetMaxima(maxima);
				minuit.SetMinima(minima);
				minuit.SetInitialValues(InitialValues);
				minuit.SetInitialErrors(InitialErrors);



				// ////////////     Now perform the fits ////////////         
				FitResult result_minuit = minuit.Optimise(&lhFunction_gSearch);

				std::vector<double> fit_minuit = result_minuit.GetBestFit();
				std::cout<<"Minuit Results, pdf 0 = bkg, 1 = sig"<<std::endl;
				result_minuit.Print();
				result_minuit.SaveAs("simpleFit_result_minuit.txt");


				if(QmetHast){
								//----------------------MetHast--------------------------
								MetropolisHastings metHast; 
								// metHast.SetInitialTrial(fit_minuit);
								metHast.SetMaxIter(100000); 
								metHast.SetMaxima(maxima); 
								metHast.SetMinima(minima);         
								metHast.SetFlipSign(true); 
								metHast.SetTestStatLogged(true); 

								std::vector<double> sigmas(inputFiles.size(),10); 
								metHast.SetSigmas(sigmas); 
								FitResult result_metHast = metHast.Optimise(&lhFunction_gSearch); 

								std::vector<double> fit_metHast =result_metHast.GetBestFit();         
								std::cout<<"MetHast Results, pdf 0 = bkg, 1 = sig"<<std::endl; 
								// std::vector<double> fit_metHast(2,1); 
								std::vector<size_t> indices; 
								indices.push_back(0);
								indices.push_back(1);
								Histogram hist = result_metHast.GetStatSpace(); 
								PdfConverter::ToTH2D(hist.Marginalise(indices)).SaveAs("lh_2d.root"); 
								// PdfConverter::ToTH1D(hist.Marginalise(0)).SaveAs("pdf_norm0.root"); 
								// PdfConverter::ToTH1D(hist.Marginalise(1)).SaveAs("pdf_norm1.root"); 
								// PdfConverter::ToTH1D(hist.Marginalise(2)).SaveAs("pdf_norm2.root"); 
								// PdfConverter::ToTH1D(hist.Marginalise(3)).SaveAs("pdf_norm3.root"); 
								result_metHast.Print(); 
								result_metHast.SaveAs("simpleFit_result_metHast.txt"); 
				}
				// ///////////////Setting up the histograms ///////////////////// 

				gStyle->SetOptStat(kFALSE); 
				std::vector<TH1D*> clone_histList;
				std::cout << "hist List length = "<<histList.size() << std::endl;

				for(int i=0;i<histList.size();i++){
								clone_histList.push_back((TH1D*) histList[i].Clone());
								std::cout << "fit_scale = "<<fit_minuit[i] << std::endl;

								if(i!=0){ clone_histList[i]->Scale(fit_minuit[i-1]);}

				}

				TH1D * complete_blank= new TH1D("","",60,0,3); 
				TH1D * complete_minuit= new TH1D("","",60,0,3); 
				TH1D * complete_metHast= new TH1D("","",60,0,3); 

				for (int i = 1; i <clone_histList.size(); i++) {
								complete_minuit->Add(clone_histList[i],1); 
				}


				complete_minuit->SetLineColor(kRed); 
				complete_minuit->SetLineWidth(2); 
				complete_minuit->SetFillColorAlpha(kRed,0.5); 
				if(QmetHast){
								for (int i = 1; i <clone_histList.size(); i++) {
												complete_metHast->Add(clone_histList[i],1); 
								}

								complete_metHast->SetLineColor(kBlack); 
								complete_metHast->SetLineWidth(2); 
								complete_metHast->SetFillColorAlpha(kBlack,0.5); 
				}


				TCanvas * fitCan= new TCanvas(); 
				complete_blank->SetTitle("Complete spectra for different methods"); 
				complete_blank->GetXaxis()->SetTitle("Energy (Mev)"); 
				complete_blank->GetYaxis()->SetTitle("Counts/ 50 keV bin"); 
				complete_blank->SetMaximum(500);
				complete_blank->Draw(); 
				complete_minuit->Draw("same"); 
				if(QmetHast) complete_metHast->Draw("same"); 

				TLegend* leg =new TLegend(0.7,0.7,0.9,0.9); 
				leg->AddEntry(complete_minuit,"Minuit","lf"); 
				if(QmetHast) leg->AddEntry(complete_metHast,"Metropolis Hastings","lf"); 
				leg->Draw(); 

				fitCan->Print("All_Fits.png"); 


				//_+_+_+_+_+_+_+_+_+_Double Panel Plots_+_+_+_+_+_+_+_+_+_ 

				// ------- Create a canvas: 
				TCanvas* diff = new TCanvas("diff","",800,800);
				diff->cd(); 
				// -------------- Top panel 
				gStyle->SetOptStat(kFALSE);  
				TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); 
				pad1->Draw(); 
				pad1->cd(); 
				pad1->SetGrid(kTRUE); 
				pad1->SetBottomMargin(0.00); 
				gPad->RedrawAxis();  

				histList[0].SetTitle("Spectrum fit with fractional bin error");
				histList[0].GetYaxis()->SetTitle("Counts / 50 keV bin");
				histList[0].GetYaxis()->SetTitleOffset(1); 
				histList[0].SetFillColorAlpha(kGreen,0.5); 
				histList[0].SetMaximum(500); 
				histList[0].Draw(); 
				complete_minuit->Draw("same e"); 
				if(QmetHast) complete_metHast->Draw("same e"); 
				leg->Draw(); 
				diff->cd(); 
				// -------------- Bottom panel 
				TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3); 
				pad2->SetTopMargin(0.00); 
				pad2->Draw(); 
				pad2->cd(); 
				pad2->SetBottomMargin(0.3); 
				pad2->SetGrid(kTRUE); 
				gStyle->SetOptStat(kFALSE);  


				// TH1D * diff_gSearch= diffHist(&comPlot,complete_gSearch); 
				TH1D * diff_minuit= diffHist(&histList[0],complete_minuit); 
				diff_minuit->SetLineColor(kRed); 
				diff_minuit->SetLineWidth(2); 
				diff_minuit->SetFillColorAlpha(kRed,0.1); 



				diff_minuit->GetXaxis()->SetTitle("Energy (MeV)"); 
				diff_minuit->GetYaxis()->SetTitle("Fractional bin error"); 
				diff_minuit->GetYaxis()->SetTitleOffset(0.5); 
				diff_minuit->GetXaxis()->SetLabelSize(0.1); 
				diff_minuit->GetXaxis()->SetTitleSize(0.1); 
				diff_minuit->GetYaxis()->SetLabelSize(0.1); 
				diff_minuit->GetYaxis()->SetTitleSize(0.1); 
				
				diff_minuit->Draw(); 

				if(QmetHast){
							 	TH1D * diff_metHast= diffHist(&histList[0],complete_metHast); 
								diff_metHast->SetLineColor(kBlack); 
								diff_metHast->SetLineWidth(2); 
								diff_metHast->SetFillColorAlpha(kBlack,0.1); 
								diff_metHast->Draw("same");
				}
				// gStyle->SetOptStat(0); 

				diff->Print("doublePlanel.png"); 
				std::cout<<" Number diff = "<< histList[0].Integral() - fit_minuit[0]-fit_minuit[1]-fit_minuit[2] <<std::endl; 

				return 0; 
}
