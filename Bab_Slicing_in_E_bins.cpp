// A simple fit in energy for signal and a background        
#include "util.h"
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <PdfConverter.h>        
#include <GridSearch.h>        
#include <TFile.h>        
#include <TF1.h>        
#include <TH1D.h>        
#include <TCanvas.h>        
#include <BoxCut.h>        
#include <BoolCut.h>        
#include <CutCollection.h>        
#include <BinnedPdfShrinker.h>        

#include <vector>
#include <cmath>

const std::string bgMCfile    = "";
const std::string sigMCfile   = "";
const std::string bgTreeName  = "";
const std::string sigTreeName = "";

const std::string dataFile = "";
const std::string dataTreeName = "";


void fillPdfs(const std::vector<std::string> betaFileList, const std::vector<std::string> alphaFileList, BinnedPdf& betaPdf, BinnedPdf& alphaPdf){


				const	std::string treeName="output";

				for (int i = 0; i < betaFileList.size(); i++) {
								std::cout<<"Filename : "<<betaFileList[i]<<std::endl;

								ROOTNtuple betaFile( betaFileList[i] , treeName );

								for(size_t i = 0; i < betaFile.GetNEntries(); i++){        
												betaPdf.Fill(betaFile.GetEntry(i));        
								}        
				} 

				for (int i = 0; i < alphaFileList.size(); i++) {
								std::cout<<"Filename : "<<alphaFileList[i]<<std::endl; 

								ROOTNtuple alphaFile( alphaFileList[i] , treeName );

								for(size_t i = 0; i < alphaFile.GetNEntries(); i++){        
												alphaPdf.Fill(alphaFile.GetEntry(i));        
								}        
				}

				std::cout << "Filled pdfs " << std::endl;        
}

void fillPdfs(const std::vector<std::string> betaFileList, const std::vector<std::string> alphaFileList, std::vector<BinnedPdf>& betaPdfs, std::vector<BinnedPdf>& alphaPdfs){

				double evIndex_cut=0;
				double fitValid_cut =1;
				double Rhigh=5000;


				const	std::string treeName="output";

				for (int i = 0; i < betaFileList.size(); i++) {
								std::cout<<"Filename : "<<betaFileList[i]<<std::endl;

								ROOTNtuple betaFile( betaFileList[i] , treeName );

								double sliceWidth=0.1;
								int counter=0;
								//It is to 2.5 because we count from 0;
								for (double lowerCut = 0; lowerCut < 2.5; lowerCut+=sliceWidth) {
												// std::cout << "lower cut = "<<lowerCut << std::endl;

												
												CutCollection cutCol;
												BoxCut Ecut(0,lowerCut,lowerCut+sliceWidth);
												cutCol.AddCut(Ecut);
												BoxCut Rcut(1,0,Rhigh);
												cutCol.AddCut(Rcut);
												BoxCut Babcut(2,-10000,10000);
												cutCol.AddCut(Babcut);
												BoolCut fitValidcut(3,fitValid_cut);
												cutCol.AddCut(fitValidcut);
												BoolCut evIndexcut(4,evIndex_cut);
												cutCol.AddCut(evIndexcut);

												for(size_t a = 0; a < betaFile.GetNEntries(); a++){        
																if (cutCol.PassesCuts(betaFile.GetEntry(a))) {
																				betaPdfs[counter].Fill(betaFile.GetEntry(a));        
																}
												}        
												counter++;
								}
				} 

				for (int i = 0; i < alphaFileList.size(); i++) {
								std::cout<<"Filename : "<<alphaFileList[i]<<std::endl; 

								ROOTNtuple alphaFile( alphaFileList[i] , treeName );

								double sliceWidth=0.1;
								int counter=0;

								for (double lowerCut = 0; lowerCut < 2.5; lowerCut+=sliceWidth) {
												// std::cout << "lower cut = "<<lowerCut << std::endl;

												CutCollection cutCol;
												BoxCut Ecut(0,lowerCut,lowerCut+sliceWidth);
												cutCol.AddCut(Ecut);
												BoxCut Rcut(1,0,Rhigh);
												cutCol.AddCut(Rcut);
												BoxCut Babcut(2,-10000,10000);
												cutCol.AddCut(Babcut);
												BoolCut fitValidcut(3,fitValid_cut);
												cutCol.AddCut(fitValidcut);
												BoolCut evIndexcut(4,evIndex_cut);
												cutCol.AddCut(evIndexcut);

												for(size_t a = 0; a < alphaFile.GetNEntries(); a++){        

																if (cutCol.PassesCuts(alphaFile.GetEntry(a))) {
																				alphaPdfs[counter].Fill(alphaFile.GetEntry(a));        
																}
												}        
												counter++;
								}

				}

				std::cout << "Filled pdfs " << std::endl;        
}


int main(){        
				////////////////////        
				// 1. Set Up PDFs //        
				////////////////////        

				// Set up binning        
				AxisCollection axes;        
				
				double Emin=0;				
				double Emax=2.6;				
				double Ebins=26;				

				double Rmin=0;				
				double Rmax=6000;				
				double Rbins=6;				

				double Babmin=-200;				
				double Babmax=200;				
				double Babbins=Babmax-Babmin;

				axes.AddAxis(PdfAxis("mcEdepQuenched", Emin, Emax, Ebins, "MCEdepQuenched"));
				axes.AddAxis(PdfAxis("mcPosr", Rmin,Rmax, Rbins, "MCPosr"));
				axes.AddAxis(PdfAxis("berkeleyAlphaBeta",Babmin,Babmax, Babbins, "BerekleyAlphaBeta"));
				axes.AddAxis(PdfAxis("fitValid",0,5,5,"FitValid"));
				axes.AddAxis(PdfAxis("evIndex",0,5,5, "EvIndex"));

				// Only interested in first bit of data ntuple        

				std::vector<size_t> indices;
				for (int i = 0; i < axes.GetNDimensions(); i++) {
								indices.push_back(i);
				}

				DataRepresentation dataRep(indices);        

				// Set up pdf with these bins in this observable        
				BinnedPdf betaPdf(axes);      betaPdf.SetDataRep(dataRep);        
				BinnedPdf alphaPdf(axes);  		alphaPdf.SetDataRep(dataRep);        

				std::vector<BinnedPdf> betaPdfs;
				std::vector<BinnedPdf> alphaPdfs;
				//This Ebins Pdfs are created because you are slicing in energy. 
				for (int i = 0; i < Ebins; ++i){

							BinnedPdf bPdf(axes);
							betaPdfs.push_back(bPdf);	
							betaPdfs[i].SetDataRep(dataRep);

							BinnedPdf aPdf(axes);
							alphaPdfs.push_back(aPdf);	
						 	alphaPdfs[i].SetDataRep(dataRep);
				}

				std::cout << "Initialised Pdfs" << std::endl;        

				std::vector<std::string> betaFileList = UTIL::glob("testData/solar/beta","electron");
				std::vector<std::string> alphaFileList = UTIL::glob("testData/solar/alpha","alpha");

				// fillPdfs(betaFileList, alphaFileList, betaPdf, alphaPdf);
				fillPdfs(betaFileList, alphaFileList, betaPdfs, alphaPdfs);


				TFile pdf_Bab_File("pdfs_Bab.root","recreate");
				pdf_Bab_File.cd();

				for (int i = 0; i < Ebins; ++i){

								TH1D Bab_beta_1 = PdfConverter::ToTH1D(betaPdfs[i].Marginalise(2),false);
								TH1D Bab_alpha_1 = PdfConverter::ToTH1D(alphaPdfs[i].Marginalise(2),false);

								Bab_beta_1.SetName(Form("Bab_beta_%d",i));
								Bab_alpha_1.SetName(Form("Bab_alpha_%d",i));

								TCanvas * c1 = new TCanvas();

								Bab_beta_1.SetMaximum(14000);
								Bab_beta_1.Draw();
								Bab_alpha_1.Draw("same");

								c1->Print(Form("Bab_pdf_%d.png",i));

								Bab_beta_1.Write();
								Bab_alpha_1.Write();
				}


				std::vector<double> norms_alpha,means_alpha, sigmas_alpha;
				std::vector<double> norms_beta,means_beta, sigmas_beta;
				for (int i = 0; i < Ebins; ++i){

								betaPdfs[i].Normalise();
								alphaPdfs[i].Normalise();

								TH1D Bab_beta_1 = PdfConverter::ToTH1D(betaPdfs[i].Marginalise(2),false);
								TH1D Bab_alpha_1 = PdfConverter::ToTH1D(alphaPdfs[i].Marginalise(2),false);

 								TF1 *fit_beta = new TF1("fit_beta","[0]*exp(-0.5*((x-[1])/[2])^2)", -200, 200);
 								TF1 *fit_alpha = new TF1("fit_alpha","[0]*exp(-0.5*((x-[1])/[2])^2)", -200, 200);

								fit_beta->SetParName(0,"Norm");
								fit_beta->SetParName(1,"Mean");
								fit_beta->SetParName(2,"Sigma");
								fit_beta->SetParameter("Norm",100);
								fit_beta->SetParameter("Mean",0);
								fit_beta->SetParameter("Sigma",1);

								fit_alpha->SetParName(0,"Norm");
								fit_alpha->SetParName(1,"Mean");
								fit_alpha->SetParName(2,"Sigma");
								fit_alpha->SetParameter("Norm",100);
								fit_alpha->SetParameter("Mean",0);
								fit_alpha->SetParameter("Sigma",1);

								Bab_beta_1.Fit("fit_beta");
								Bab_alpha_1.Fit("fit_alpha");

								norms_beta.push_back(fit_beta->GetParameter(0));
								means_beta.push_back(fit_beta->GetParameter(1));
								sigmas_beta.push_back(fit_beta->GetParameter(2));

								norms_alpha.push_back(fit_alpha->GetParameter(0));
								means_alpha.push_back(fit_alpha->GetParameter(1));
								sigmas_alpha.push_back(fit_alpha->GetParameter(2));

								Bab_beta_1.SetName(Form("Bab_beta_%d_WithFit",i));
								Bab_alpha_1.SetName(Form("Bab_alpha_%d_WithFit",i));

								TCanvas * c1 = new TCanvas();

								Bab_beta_1.SetMaximum(14000);
								Bab_beta_1.Draw();
								Bab_alpha_1.Draw("same");

								c1->Print(Form("Bab_pdf_%d_with_fit.png",i));

								Bab_beta_1.Write();
								Bab_alpha_1.Write();

				}

				

				pdf_Bab_File.Close();




				return 0;


				/////////////////////////////////////        
				// 2. Fill with data and normalise //        
				/////////////////////////////////////        

				BinnedPdf bgPdf(axes);      bgPdf.SetDataRep(dataRep);        
				BinnedPdf signalPdf(axes);  signalPdf.SetDataRep(dataRep);        
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



				////////////////////////////        
				// 3. Set Up LH function  //        
				////////////////////////////        
				ROOTNtuple dataNt(dataFile, dataTreeName);
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
