#include <catch.hpp>
#include <BinnedNLLH.h>
#include <Gaussian.h>
#include <PdfConverter.h>
#include <OXSXDataSet.h>
#include <iostream>

TEST_CASE("Binned NLLH, 3 rates no systematics"){
    Gaussian gaus1(1, 10);
    Gaussian gaus2(1, 20);
    Gaussian gaus3(1, 30);

    AxisCollection axes;
    axes.AddAxis(PdfAxis("axis1", -40, 40 , 200));

    BinnedPdf pdf1 = PdfConverter::ToBinnedPdf(gaus1, axes);
    BinnedPdf pdf2 = PdfConverter::ToBinnedPdf(gaus2, axes);
    BinnedPdf pdf3 = PdfConverter::ToBinnedPdf(gaus3, axes);


    size_t centralBin = pdf1.FindBin(std::vector<double>(1,0));
    double prob1 = pdf1.GetBinContent(centralBin);
    double prob2 = pdf2.GetBinContent(centralBin);
    double prob3 = pdf3.GetBinContent(centralBin);

    pdf1.SetDataRep(0);
    pdf2.SetDataRep(0);
    pdf3.SetDataRep(0);
    
    BinnedNLLH lh;
    lh.AddPdf(pdf1);
    lh.AddPdf(pdf2);
    lh.AddPdf(pdf3);

    OXSXDataSet data;
    data.AddEntry(EventData(std::vector<double>(1, 0)));
    
    lh.SetDataSet(&data);
    
    lh.RegisterFitComponents();
    SECTION("Correct Probability"){
        double sumLogProb = -log(prob1 + prob2 + prob3);
        double sumNorm    = 3;
        
        lh.SetParameters(std::vector<double>(3, 1));
        REQUIRE(lh.Evaluate() == Approx(sumNorm + sumLogProb));
    }

    SECTION("Correct Probability with constraint"){
        lh.SetConstraint("Pdf Normalisation 0", 3, 1);

        double sumLogProb = -log(prob1 + prob2 + prob3);
        double sumNorm    = 3;
        double constraint = 2;

        lh.SetParameters(std::vector<double>(3, 1));
        REQUIRE(lh.Evaluate() == Approx(sumNorm + sumLogProb + constraint));
    }
}
