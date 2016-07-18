// A simple fit in energy for signal and a background        
#include "util.h"
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <PdfConverter.h>        
#include <GridSearch.h>        
#include <TFile.h>        
#include <TH1D.h>        
#include <TCanvas.h>        
#include <BoxCut.h>        
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


				const	std::string treeName="output";

				for (int i = 0; i < betaFileList.size(); i++) {
								std::cout<<"Filename : "<<betaFileList[i]<<std::endl;

								ROOTNtuple betaFile( betaFileList[i] , treeName );

								double sliceWidth=0.1;
								int counter=0;
								for (double lowerCut = 0; lowerCut < 2.6; lowerCut+=sliceWidth) {
												// std::cout << "lower cut = "<<lowerCut << std::endl;

												
												CutCollection cutCol;
												BoxCut Ecut(0,lowerCut,lowerCut+sliceWidth);
												cutCol.AddCut(Ecut);
												BoxCut Rcut(1,0,10000);
												cutCol.AddCut(Rcut);
												BoxCut Babcut(2,-10000,10000);
												cutCol.AddCut(Babcut);

												// std::cout << "betaFile entries = "<<  betaFile.GetNEntries()<< std::endl;
												for(size_t a = 0; a < betaFile.GetNEntries(); a++){        
																// betaPdf.Fill(betaFile.GetEntry(a));        
																if (cutCol.PassesCuts(betaFile.GetEntry(a))) {
																				// std::cout << "betaFile inside the loop = "<< std::endl;
																				betaPdfs[counter].Fill(betaFile.GetEntry(a));        
																}
																// betaPdfs[counter].Fill(betaFile.GetEntry(i));        
																// std::cout << "counter = "<< counter << std::endl;
												}        
												counter++;
								}
				} 

				for (int i = 0; i < alphaFileList.size(); i++) {
								std::cout<<"Filename : "<<alphaFileList[i]<<std::endl; 

								ROOTNtuple alphaFile( alphaFileList[i] , treeName );

								double sliceWidth=0.1;
								int counter=0;

								for (double lowerCut = 0; lowerCut < 2.6; lowerCut+=sliceWidth) {
												// std::cout << "lower cut = "<<lowerCut << std::endl;

												CutCollection cutCol;
												BoxCut Ecut(0,lowerCut,lowerCut+sliceWidth);
												cutCol.AddCut(Ecut);
												BoxCut Rcut(1,0,10000);
												cutCol.AddCut(Rcut);
												BoxCut Babcut(2,-10000,10000);
												cutCol.AddCut(Babcut);

												// std::cout << "alphaFile entries = "<<  alphaFile.GetNEntries()<< std::endl;
												for(size_t a = 0; a < alphaFile.GetNEntries(); a++){        

																if (cutCol.PassesCuts(alphaFile.GetEntry(a))) {
																				// std::cout << "alphaFile inside loop"<< std::endl;
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

				std::cout << "beta PDFs 1 = "<<betaPdfs[1].Means()[1] << std::endl;
				std::cout << "alpha PDFs 1 = "<<alphaPdfs[1].Means()[1] << std::endl;

				for (int i = 0; i < Ebins; ++i){

								TH1D Bab_beta_1 = PdfConverter::ToTH1D(betaPdfs[i].Marginalise(2),false);
								TH1D Bab_alpha_1 = PdfConverter::ToTH1D(alphaPdfs[i].Marginalise(2),false);

								Bab_beta_1.SetName(Form("Bab_beta_%d",i));
								Bab_alpha_1.SetName(Form("Bab_alpha_%d",i));

								TCanvas * c1 = new TCanvas();

								Bab_beta_1.Draw();
								Bab_alpha_1.Draw("same");

								// c1->Print("Bab_pdf_1.png");
								c1->Print(Form("Bab_pdf_%d.png",i));

								Bab_beta_1.Write();
								Bab_alpha_1.Write();
				}

				pdf_Bab_File.Close();
				return 0;
				// betaPdf.Normalise();        
				// alphaPdf.Normalise();        

				TH1D betaPdfHist_energy = PdfConverter::ToTH1D(betaPdf.Marginalise(0),false);
				TH1D alphaPdfHist_energy = PdfConverter::ToTH1D(alphaPdf.Marginalise(0),false);
				TH1D betaPdfHist_Bab = PdfConverter::ToTH1D(betaPdf.Marginalise(2),false);
				TH1D alphaPdfHist_Bab = PdfConverter::ToTH1D(alphaPdf.Marginalise(2),false);
				betaPdfHist_energy.SetName("mcEdepQuenched_Beta");
				alphaPdfHist_energy.SetName("mcEdepQuenched_Alpha");
				betaPdfHist_Bab.SetName("Bab_Beta");
				alphaPdfHist_Bab.SetName("Bab_Alpha");
				TFile pdfFile("pdfs.root","update");
				pdfFile.cd();
				betaPdfHist_energy.Write();
				alphaPdfHist_energy.Write();
				betaPdfHist_Bab.Write();
				alphaPdfHist_Bab.Write();
				pdfFile.Close();


				
				BinnedPdfShrinker shrinker;
				shrinker.SetBuffer(0,10,10);
				BinnedPdf mmm = shrinker.ShrinkPdf(betaPdf);

				TH1D ASkrink= PdfConverter::ToTH1D(mmm.Marginalise(0),false);
				TH1D ASkrink_2= PdfConverter::ToTH1D(mmm.Marginalise(2),false);

				ASkrink.SetName("Energy After skrink");
				ASkrink_2.SetName("Bab After skrink");

				TFile shrinkFile("pdfs_after_Shrink.root","update");
				shrinkFile.cd();
				ASkrink.Write();
				ASkrink_2.Write();
				shrinkFile.Close();
				

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
