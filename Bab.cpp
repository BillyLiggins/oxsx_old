// A simple fit in energy for signal and a background        
#include "util.h"
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <PdfConverter.h>        
#include <GridSearch.h>        
#include <TFile.h>        
#include <TH1D.h>        

#include <vector>

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


int main(){        
				////////////////////        
				// 1. Set Up PDFs //        
				////////////////////        

				// Set up binning        
				AxisCollection axes;        
				axes.AddAxis(PdfAxis("mcEdepQuenched", 0, 2.6, 26, "MCEdepQuenched"));
				axes.AddAxis(PdfAxis("mcPosr", 0, 6000, 6, "MCPosr"));
				axes.AddAxis(PdfAxis("berkeleyAlphaBeta",-200, 200, 400, "BerekleyAlphaBeta"));

				// Only interested in first bit of data ntuple        

				std::vector<size_t> indices;
				indices.push_back(0);
				indices.push_back(1);
				indices.push_back(2);
				DataRepresentation dataRep(indices);        

				// Set up pdf with these bins in this observable        
				BinnedPdf betaPdf(axes);      betaPdf.SetDataRep(dataRep);        
				BinnedPdf alphaPdf(axes);  alphaPdf.SetDataRep(dataRep);        

				std::cout << "Initialised Pdfs" << std::endl;        

				std::vector<std::string> betaFileList = UTIL::glob("testData/solar/beta","electron");
				std::vector<std::string> alphaFileList = UTIL::glob("testData/solar/alpha","alpha");

				fillPdfs(betaFileList, alphaFileList, betaPdf, alphaPdf);


				betaPdf.Normalise();        
				alphaPdf.Normalise();        

				TH1D betaPdfHist_energy = PdfConverter::ToTH1D(betaPdf.Marginalise(0),false);
				TH1D alphaPdfHist_energy = PdfConverter::ToTH1D(alphaPdf.Marginalise(0),false);
				betaPdfHist_energy.SetName("mcEdepQuenched_Beta");
				alphaPdfHist_energy.SetName("mcEdepQuenched_Alpha");
				TFile pdfFile("pdfs.root","update");
				pdfFile.cd();
				betaPdfHist_energy.Write();
				alphaPdfHist_energy.Write();
				pdfFile.Close();



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
