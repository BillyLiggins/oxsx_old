// A simple fit in energy for signal and a background        
#include <BinnedPdf.h>        
#include <ROOTNtuple.h>        
#include <BinnedNLLH.h>        
#include <GridSearch.h>        
#include <DataSetGenerator.h>        
#include <OXSXDataSet.h>         
#include <PdfConverter.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <Rand.h>
#include <Gaussian.h>
#include <Convolution.h>
#include <Scale.h>
#include <Minuit.h>
#include <IntegrableFunction.h> 
#include <BoxCut.h>

// const std::string sig1MCfile   = "testData/Bi210/SolarBi210_complete.ntuple_oxsx.root";
// const std::string sig2MCfile   = "testData/Po210/SolarPo210_complete.ntuple_oxsx.root";
const std::string sig1MCfile   = "testData/Bi210/SolarBi210_data.ntuple_oxsx.root";
const std::string sig2MCfile   = "testData/Po210/SolarPo210_data.ntuple_oxsx.root";
const std::string sig1TreeName = "output";
const std::string sig2TreeName = "output";

const std::string dataFile1 = "testData/Bi210/SolarBi210_data.ntuple_oxsx.root";
const std::string dataFile2 = "testData/Po210/SolarPo210_data.ntuple_oxsx.root";
const std::string dataTreeName = "output";

int main(){

		bool addsys=true;
		// start by setting the random number generator - this makes it change each time run.
		Rand::SetSeed(0);

		////////////////////        
		// 1. Set Up PDFs //        
		////////////////////        

		// Set up binning        
		AxisCollection axes;        
		axes.AddAxis(PdfAxis("nhits", 0,600,30, "nhits"));

		// Only interested in first bit of data ntuple        
		DataRepresentation dataRep(0);        

		// Set up pdf with these bins in this observable        
		BinnedPdf sig1Pdf(axes);      sig1Pdf.SetDataRep(dataRep);
		BinnedPdf sig2Pdf(axes);      sig2Pdf.SetDataRep(dataRep);

		std::cout << "Initialised Pdfs" << std::endl;        

		/////////////////////////////////////        
		// 2. Fill with data and normalise //        
		/////////////////////////////////////        

		ROOTNtuple sig1MC(sig1MCfile, sig1TreeName);
		ROOTNtuple sig2MC(sig2MCfile, sig2TreeName);

		BoxCut boxCut(0,50,550); // cuts away the overflows

		for(size_t i = 0; i < sig1MC.GetNEntries(); i++){
				if(boxCut.PassesCut(sig1MC.GetEntry(i))){
						sig1Pdf.Fill(sig1MC.GetEntry(i));
				}
		}        

		for(size_t i = 0; i < sig2MC.GetNEntries(); i++){
				if(boxCut.PassesCut(sig2MC.GetEntry(i))){
						sig2Pdf.Fill(sig2MC.GetEntry(i));
				}
		}        

		sig1Pdf.Normalise();
		sig2Pdf.Normalise();

		std::cout << "Filled pdfs " << std::endl;

		// lets have a look at them
		TH1D th1fromPdf   = PdfConverter::ToTH1D(sig1Pdf, false);
		th1fromPdf.SetName("pdf_1");
		TH1D th2fromPdf   = PdfConverter::ToTH1D(sig2Pdf, false);
		th2fromPdf.SetName("pdf_2");
		TFile CheckPDFs("CheckPdfs.root","recreate");
		th1fromPdf.Write();
		th2fromPdf.Write();

		// Systematics...
		// Now smear with a gaussian nhits resolution
		Gaussian gaussianPdf(0,20);   // Set mean and sigma
		//    gaussianPdf.MakeFittable();   // redundant - means could set parameters through a set param command
		Convolution smearer;
		smearer.SetFunction(&gaussianPdf);
		smearer.SetAxes(axes);        // axes as data - just 1 in nhits
		smearer.SetDataRep(dataRep);
		smearer.SetPdfDataRep(dataRep);
		smearer.MakeFittable();
		smearer.Construct();

		// check the smearer is getting the right parameters
		std::cout << smearer.GetParameters()[0] << " " << smearer.GetParameters()[1]<< std::endl;

		// Here I want to take a look at the PDFs smeared
		BinnedPdf sig1PdfSmear = smearer(sig1Pdf);
		BinnedPdf sig2PdfSmear = smearer(sig2Pdf);

		TH1D th1fromPdf_smear   = PdfConverter::ToTH1D(sig1PdfSmear, false);
		th1fromPdf_smear.SetName("pdf_1_smear");
		TH1D th2fromPdf_smear   = PdfConverter::ToTH1D(sig2PdfSmear, false);
		th2fromPdf_smear.SetName("pdf_2_smear");
		TCanvas * normCan = new TCanvas();
		th1fromPdf_smear.Draw();
		th2fromPdf_smear.Draw("same");
		normCan->Print("Test.png");
		th1fromPdf_smear.Write();
		th2fromPdf_smear.Write();
		CheckPDFs.Close();

		// Scale scaler ;
		// scaler.SetAxes(axes);        // axes as data - just 1 in nhits
		// scaler.SetDataRep(dataRep);
		// scaler.SetPdfDataRep(dataRep);
		// scaler.SetScaleFactor(1.);
		// scaler.MakeFittable();
		// scaler.Construct();

		////////////////////////////        
		// 3. Set Up LH function  //        
		////////////////////////////        
		ROOTNtuple bi210(dataFile1, dataTreeName);
		ROOTNtuple po210(dataFile2, dataTreeName);

		// Set up a loop to do bias tests
		//    const int ntests = 1;
		const int ntests = 100;
		int nA = 300;
		int nB = 500;

		CutCollection thecuts;
		thecuts.AddCut(boxCut);

		std::vector<double> fits[ntests];
		for(int itest=0;itest<ntests;++itest){

				// DataSetGenerator
				DataSetGenerator dataGen;
				dataGen.AddDataSet(&bi210, nA);  // second arg is the expected rate
				dataGen.AddDataSet(&po210, nB);
				dataGen.SetCuts(thecuts);

				// Generate a data set
				OXSXDataSet fakeData  = dataGen.PoissonFluctuatedDataSet();

				// Try looking at this data set
				TFile CheckPDFs2("CheckPdfs.root","UPDATE");
				TH1D thedata(Form("thedata_%d",itest),Form("Data histogram for set %d", itest),24,0,600);
				for(int id=0;id<fakeData.GetNEntries();++id){
						float val = fakeData.GetEntry(id).GetDatum(0);
						thedata.Fill(val);
				}
				thedata.Write();

				//    OXSXDataSet fakeData  = dataGen.ExpectedRatesDataSet();

				// Set up the binned fit for NLLH
				/*    BinnedPdfManager pManager;
					  pManager.AddPdf(sig1Pdf);
					  pManager.AddPdf(sig2Pdf);

					  SystematicManager sManager;
					  sManager.Add(&smearer);

					  BinnedNLLH lh;
					  lh.SetPdfManager(pManager);
					  lh.SetSystematicManager(sManager);
					  lh.SetDataSet(&fakeData);
					  */
				BinnedNLLH lh;
				lh.SetBufferAsOverflow(false); // when the buffer is applied it throws away the bins that are smeared off the end of the buffer and renormalises, rather than piling the probability into the first bin.
				lh.SetBuffer(0,2,2);    // creates a buff of 2 bins at either side of axis 0
				lh.AddCut(boxCut);      // need to make sure data is binned with this cut as well as MC
				lh.SetDataSet(&fakeData);
				lh.AddPdf(sig1Pdf);
				lh.AddPdf(sig2Pdf);
				if(addsys)lh.AddSystematic(&smearer);//lh.AddSystematic(&scaler);

				std::cout << "Built LH function for test " << itest << std::endl;

				// Set up the optimisation
				std::vector<double> minima;
				minima.push_back(10);
				minima.push_back(10);
				if(addsys)  minima.push_back(-10.); // Gaussian mean and sigma
				if(addsys)  minima.push_back(0.);
				// if(addsys)  minima.push_back(1.); //Scale 
				std::vector<double> maxima;
				maxima.push_back(1000);
				maxima.push_back(1000);
				if(addsys)  maxima.push_back(10.);
				if(addsys)  maxima.push_back(2.);
				// if(addsys)  maxima.push_back(4.);
				std::vector<double> initialval;
				initialval.push_back(200);
				initialval.push_back(200);
				if(addsys)  initialval.push_back(-2.);
				if(addsys)  initialval.push_back(5.);
				// if(addsys)  initialval.push_back(2.5);
				std::vector<double> initialerr;
				initialerr.push_back(20);
				initialerr.push_back(20);
				if(addsys)  initialerr.push_back(0.1);
				if(addsys)  initialerr.push_back(0.1);
				// if(addsys)  initialerr.push_back(0.1);

				Minuit min;
				min.SetMinima(minima);
				min.SetMaxima(maxima);
				min.SetInitialValues(initialval);
				min.SetInitialErrors(initialerr);

				////////////
				// 4. Fit //
				////////////
				std::cout << "Starting Fit for test " << itest << std::endl;

				FitResult result = min.Optimise(&lh);

				fits[itest] = result.GetBestFit();
				result.Print();
		}

		// Now do something to plot all the fits.
		TH1D *hfitsA = new TH1D("hfitsA","",100,0,nA*1.2);
		TH1D *hfitsB = new TH1D("hfitsB","",100,nB*0.8,nB*1.2);
		TH1D *hfitsC = new TH1D("hfitsC","",100,-10,10);
		TH1D *hfitsD = new TH1D("hfitsD","",100,-10,10);
		// TH1D *hfitsE = new TH1D("hfitsE","",100,-2,2);
		float num[ntests];
		for(int i=0;i<ntests;++i){
				hfitsA->Fill(fits[i][0]);
				hfitsB->Fill(fits[i][1]);
				hfitsC->Fill(fits[i][2]);
				hfitsD->Fill(fits[i][3]);
				// hfitsE->Fill(fits[i][4]);
				std::cout << fits[i][0] << " " << fits[i][1] << " " << fits[i][2] << " " << fits[i][3] << " " << fits[i][3] <<   std::endl;
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
		// Can->cd(5);
		// hfitsE->Draw();
		Can->Print("Plots.png");
		return 0;
}
