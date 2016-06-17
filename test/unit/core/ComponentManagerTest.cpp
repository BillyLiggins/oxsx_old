#include <catch.hpp>
#include <ComponentManager.h>
#include <SpectralFitPdf.h>

TEST_CASE("Stand alone component manager"){
    // some fake objects
    AxisCollection axes;
    axes.AddAxis(PdfAxis("", 0, 10, 100));
    SpectralFitPdf pdf1(axes);
    SpectralFitPdf pdf2(axes);

    ComponentManager cmpMan;

    SECTION("initialised correctly"){
        REQUIRE(cmpMan.GetTotalParameterCount() == 0);
        REQUIRE(cmpMan.GetParameters() == std::vector<double>(0));
        REQUIRE(cmpMan.GetParameterNames() == std::vector<std::string> (0));
    }
    
    SECTION("adding a component"){
        cmpMan.AddComponent(&pdf1);
        REQUIRE(cmpMan.GetTotalParameterCount() == 100);
    }

    SECTION("adding two"){
        cmpMan.AddComponent(&pdf1);
        cmpMan.AddComponent(&pdf2);
        REQUIRE(cmpMan.GetTotalParameterCount() == 100 * 2);        

        cmpMan.SetParameters(std::vector<double> (100 * 2, 8));
        REQUIRE(pdf1.GetParameters() == std::vector<double>(100, 8));
        REQUIRE(pdf2.GetParameters() == std::vector<double>(100, 8));
    }
    SECTION("clearing"){
        cmpMan.Clear();
        REQUIRE(cmpMan.GetParameters() == std::vector<double>(0));    
        REQUIRE(cmpMan.GetTotalParameterCount() == 0);
        REQUIRE(cmpMan.GetParameterNames() == std::vector<std::string>(0));
    }
    SECTION("getting parameter by name"){
        cmpMan.AddComponent(&pdf1);
        REQUIRE(cmpMan.GetParameter("Spectral fit bin 0") == 0);
        pdf1.SetBinContent(0, 10);
        REQUIRE(cmpMan.GetParameter("Spectral fit bin 0") == 10);
        
    }
}
