// // a simple fit in energy for 2n2b, po210 and bi210 backgrounds.        
// this is very messy you should make it so that the everything is in vectors.

// #include "util.h"
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

template<typename T> void printElement(T t, const int& width)
{
				 const char separator    = ' ';
				std::cout << std::left << std::setw(width) << std::setfill(separator) << t;
}

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

int main(){        
				bool QmetHast=false;
				// bool QmetHast=true;
				bool QSys=false;
				// bool QSys=true;
				Rand::SetSeed(0);
				std::vector<std::string> inputFiles;
				inputFiles.push_back("testData/TeLoaded/Bi210/complete_Bi210.ntuple_oxsx.root");
				inputFiles.push_back("testData/TeLoaded/Po210/complete_Po210.ntuple_oxsx.root");
				inputFiles.push_back("testData/TeLoaded/2n2b/complete_2n2b.ntuple_oxsx.root");
				inputFiles.push_back("testData/TeLoaded/C14/complete_C14.ntuple_oxsx.root");
				// inputFiles.push_back("testData/TeLoaded/0n2b/complete_0n2b.ntuple_oxsx.root");

				std::vector<std::string> names;
				names.push_back("data");
				names.push_back("Bi210");
				names.push_back("Po210");
				names.push_back("2n2b");
				names.push_back("C14");
				// names.push_back("0n2b");

				std::vector<std::string> inputTrees;
				inputTrees.push_back("output");
				inputTrees.push_back("output");
				inputTrees.push_back("output");
				inputTrees.push_back("output");
				// inputTrees.push_back("output");

				////////////////////        1. Set Up PDFs //      //////////////////        

				// Set up binning        
				AxisCollection axes;
				axes.AddAxis(PdfAxis("energy", 0, 3, 60,"Energy"));

				// Only interested in first bit of data ntuple        
				DataRepresentation dataRep(0);        


				// Set up pdf with these bins in this observable        


				/******************************************************************************************
				 *************************** 2. Fill with data and normalise ******************************
				 ******************************************************************************************/        

				double daysOfData =0.001;
				std::vector<double> rates;
				rates.push_back(1000); //Bi210
				rates.push_back(1000); //Po210
				rates.push_back(1000); //2n2b
				rates.push_back(1000); //C14
				// rates.push_back(1000); //0n2b

				//rates for a year of data. From DocDb <find link>
				double Bi210Rate=427368;
				double Po210Rate=1.74e7;
				double C14Rate=3.78e9;
				double _2n2bRate=find2n2bRate();

				// std::vector<double> RATES;
				// RATES.push_back(Bi210Rate);
				// RATES.push_back(Po210Rate);
				// RATES.push_back(C14Rate);
				// RATES.push_back(_2n2bRate);
				// double highestRate=*max_element(RATES.begin(),RATES.end());
				// rates=normRates(RATES, scaleForTime(highestRate,daysOfData));
				
				// rates.push_back(scaleForTime(Bi210Rate,daysOfData)); //Bi210
				// rates.push_back(scaleForTime(Po210Rate,daysOfData)); //Po210
				// rates.push_back(scaleForTime(_2n2bRate,daysOfData)); //2n2b
				// rates.push_back(scaleForTime(C14Rate,daysOfData)); //C14

				std::cout<<"Bi "<<scaleForTime(Bi210Rate,daysOfData)<<std::endl; //Bi210
				std::cout<<"Po "<<scaleForTime(Po210Rate,daysOfData)<<std::endl; //Po210
				std::cout<<"2b2n "<<scaleForTime(_2n2bRate,daysOfData)<<std::endl; //2n2b not known
				std::cout<<"c14 "<<scaleForTime(C14Rate,daysOfData)<<std::endl; //C14

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

				int ntest=100;
				std::vector<double> fits[ntest];
				for (int test = 0; test < ntest; test++) {

								BinnedPdf DataPdf(axes);
								DataPdf.SetDataRep(dataRep);        
								OXSXDataSet fakeData= dataGen.ExpectedRatesDataSet();
								// OXSXDataSet fakeData= dataGen.PoissonFluctuatedDataSet();
								for(size_t i = 0; i < fakeData.GetNEntries(); i++){        
												//Cut between 0 and 3.
												if(boxCut.PassesCut(fakeData.GetEntry(i))){
																DataPdf.Fill(fakeData.GetEntry(i));        
												}
								}        


								TH1D dataPlot= PdfConverter::ToTH1D(DataPdf,false);
								// histList.insert(histList.begin(), PdfConverter::ToTH1D(DataPdf,false));



								//Setting up the systematics 

								Convolution conv;
								Gaussian gaus(0,1); 
								conv.SetFunction(&gaus);
								conv.SetAxes(axes);
								DataRepresentation DataRep(0);
								conv.SetDataRep(DataRep);
								conv.SetPdfDataRep(DataRep);



								///////////////////// Setting up the lhFunctions ////////////////////////

								BinnedNLLH lh;
								if (QSys){
												lh.SetBufferAsOverflow(false);
												lh.SetBuffer(0,5,5);
								}
								lh.SetDataPdf(DataPdf); // initialise with the data set
								for(int i=0;i<binnedPDFList.size();i++){
												lh.AddPdf(binnedPDFList[i]);

								}
								if (QSys) lh.AddSystematic(&conv);
								std::cout << "Built LH functions " << std::endl;


								/////////////////// Set up the optimisations ////////////////////////////
								std::vector<double> minima;
								std::vector<double> maxima;
								std::vector<double> stepsizes;
								std::vector<double> InitialValues;
								std::vector<double> InitialErrors;
								for(int i=0;i<binnedPDFList.size();i++){
												minima.push_back(0);
												maxima.push_back(12000);
												stepsizes.push_back(1);
												InitialValues.push_back(900);
												InitialErrors.push_back(1);
								}
								// Set up the optimisation
								if(QSys){
												minima.push_back(-0.002); //low Gaussian mean 
												minima.push_back(0.); 	//low gaussian sigma 

												maxima.push_back(0.002);	//high Gaussian mean 
												maxima.push_back(0.5); 	//high gaussian sigma

												InitialValues.push_back(0.);
												InitialValues.push_back(0.);

												InitialErrors.push_back(0.1);
												InitialErrors.push_back(0.1);
								}

								Minuit minuit;

								minuit.SetMaxima(maxima);
								minuit.SetMinima(minima);
								minuit.SetInitialValues(InitialValues);
								minuit.SetInitialErrors(InitialErrors);



								// ////////////     Now perform the fits ////////////         
								std::cout << "Now starting optimisation." << std::endl;
								FitResult result_minuit = minuit.Optimise(&lh);

								fits[test] = result_minuit.GetBestFit();

								std::vector<double> fit_minuit = result_minuit.GetBestFit();
								std::cout<<"Minuit Results, pdf 1 = bkg, 1 = sig"<<std::endl;
								result_minuit.Print();
								result_minuit.SaveAs("simpleFit_result_minuit.txt");

								std::vector<double> fit_metHast;
								if(QmetHast){
												//----------------------MetHast--------------------------
												MetropolisHastings metHast; 
												// metHast.SetInitialTrial(fit_minuit);
												metHast.SetMaxIter(100000); 
												metHast.SetMaxima(maxima); 
												metHast.SetMinima(minima);         
												metHast.SetFlipSign(true); 
												metHast.SetTestStatLogged(true); 

												std::vector<double> sigmas; 
												sigmas.push_back(1);
												sigmas.push_back(1);
												sigmas.push_back(1);
												sigmas.push_back(1);
												if(QSys)	sigmas.push_back(0.001);
												if(QSys)	sigmas.push_back(0.001);
												// std::vector<double> sigmas(inputFiles.size(),10); 
												metHast.SetSigmas(sigmas); 
												FitResult result_metHast = metHast.Optimise(&lh); 

												fit_metHast =result_metHast.GetBestFit();         
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
								std::vector<TH1D*> clone_histList_minuit, clone_histList_metHast;
								std::cout << "hist List length = "<<histList.size() << std::endl;

								for(int i=0;i<histList.size();i++){
												clone_histList_minuit.push_back((TH1D*) histList[i].Clone());
												clone_histList_metHast.push_back((TH1D*) histList[i].Clone());
												std::cout << "fit_scale = "<<fit_minuit[i] << std::endl;

												clone_histList_minuit[i]->Scale(fit_minuit[i]);
												if(QmetHast){ clone_histList_metHast[i]->Scale(fit_metHast[i]);}

								}

								TH1D * complete_blank= new TH1D("complete_blank","",60,0,3); 
								TH1D * complete_minuit= new TH1D("complete_minuit","",60,0,3); 
								TH1D * complete_metHast= new TH1D("complete_metHast","",60,0,3); 

								for (int i = 1; i <clone_histList_minuit.size(); i++) {
												complete_minuit->Add(clone_histList_minuit[i],1); 
								}


								complete_minuit->SetLineColor(kRed); 
								complete_minuit->SetLineWidth(2); 
								complete_minuit->SetFillColorAlpha(kRed,0.5); 
								if(QmetHast){
												for (int i = 1; i <clone_histList_metHast.size(); i++) {
																complete_metHast->Add(clone_histList_metHast[i],1); 
												}

												complete_metHast->SetLineColor(kBlack); 
												complete_metHast->SetLineWidth(2); 
												complete_metHast->SetFillColorAlpha(kBlack,0.5); 
								}


								TCanvas * fitCan= new TCanvas(); 
								complete_blank->SetTitle("Complete spectra for different methods"); 
								complete_blank->GetXaxis()->SetTitle("Energy (Mev)"); 
								complete_blank->GetYaxis()->SetTitle("Counts/ 50 keV bin"); 
								complete_blank->SetMaximum(5000);
								complete_blank->Draw(); 
								complete_minuit->Draw("same"); 
								if(QmetHast) complete_metHast->Draw("same"); 

								TLegend* leg =new TLegend(0.7,0.7,0.9,0.9); 
								leg->AddEntry(complete_minuit,"Minuit","lf"); 
								if(QmetHast) leg->AddEntry(complete_metHast,"Metropolis Hastings","lf"); 
								leg->Draw(); 



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

								dataPlot.SetTitle("Spectrum fit with fractional bin error");
								dataPlot.GetYaxis()->SetTitle("Counts / 50 keV bin");
								dataPlot.GetYaxis()->SetTitleOffset(1); 
								dataPlot.SetFillColorAlpha(kGreen,0.5); 
								dataPlot.SetMaximum(2000); 
								dataPlot.Draw(); 
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
								TH1D * diff_minuit= diffHist(&dataPlot,complete_minuit); 
								diff_minuit->SetLineColor(kRed); 
								diff_minuit->SetMarkerStyle(3); 
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
												TH1D * diff_metHast= diffHist(&dataPlot,complete_metHast); 
												diff_metHast->SetLineColor(kBlack); 
												diff_metHast->SetLineWidth(2); 
												diff_metHast->SetFillColorAlpha(kBlack,0.1); 
												diff_metHast->Draw("same");
								}
								// gStyle->SetOptStat(0); 



								TFile CheckPDFs2("doublePanels.root","UPDATE");

								diff->Write();
								dataPlot.Write(); 
								complete_minuit->Write("same e"); 
								if(QmetHast) complete_metHast->Write("same e"); 
								CheckPDFs2.Close();


				}

				double nA=1000;
				double nB=1000;
				double nC=1000;
				double nD=1000;

				TH1D *hfitsA = new TH1D("hfitsA","",500,nA*0.0,nA*1.2);
				TH1D *hfitsB = new TH1D("hfitsB","",500,nB*0.0,nB*1.2);
				TH1D *hfitsC = new TH1D("hfitsC","",500,nC*0.0,nC*1.2);
				TH1D *hfitsD = new TH1D("hfitsD","",500,nD*0.0,nD*1.2);
				float num[ntest];
				for(int i=0;i<ntest;++i){
								hfitsA->Fill(fits[i][0]);
								hfitsB->Fill(fits[i][1]);
								hfitsC->Fill(fits[i][2]);
								hfitsD->Fill(fits[i][3]);
								std::cout << fits[i][0] << " " << fits[i][1] << " " << fits[i][2] << " " << fits[i][3] << " " << std::endl;
								num[i] = i+1;
				}
				TCanvas *Can = new TCanvas("Can");
				Can->Divide(2,2);
				Can->cd(1);
				hfitsA->Draw();
				Can->cd(2);
				hfitsB->Draw();
				Can->cd(3);
				hfitsC->Draw();
				Can->cd(4);
				hfitsD->Draw();
				Can->Print("Plots.png");
				// printElement(2.2445335,1);
				return 0; 
}
