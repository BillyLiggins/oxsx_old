#include <FitResult.h>
#include <Histogram.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <Exceptions.h>
#include <Formatter.hpp>

FitResult::~FitResult(){
    delete fStatSpace;
}

void
FitResult::SetBestFit(const std::vector<double>& bestFit_){
    fBestFit = bestFit_;
}

std::vector<double>
FitResult::GetBestFit() const{
    return fBestFit;
}

void
FitResult::SetStatSpace(const Histogram& statSpace_){
    fStatSpace = new Histogram(statSpace_);
}

const Histogram&
FitResult::GetStatSpace() const{
    return *fStatSpace;
}

FitResult
FitResult::operator=(const FitResult& other_){
    if (!fStatSpace)
        fStatSpace = NULL;
    else
        fStatSpace = new Histogram(*other_.fStatSpace);
    fIsValid = other_.fIsValid;
    fBestFit = other_.fBestFit;
    fStatSample = other_.fStatSample;
	fParameterNames = other_.fParameterNames;
    return *this;
}

FitResult::FitResult(const FitResult& other_){
    if (!other_.fStatSpace)
        fStatSpace = NULL;

    else
        fStatSpace = new Histogram(*other_.fStatSpace);
	fParameterNames = other_.fParameterNames;
    fStatSample = other_.fStatSample;
    fBestFit = other_.fBestFit;
    fIsValid = other_.fIsValid;
}

const std::vector<std::vector<double> >&
FitResult::GetStatSample() const{
    return fStatSample;
}

void
FitResult::SetStatSample(const std::vector<std::vector<double> >& statSample_){
    fStatSample = statSample_;
}

bool
FitResult::GetValid() const {
    return fIsValid;
}

void
FitResult::SetValid(bool b_){
    fIsValid = b_;
}

void
FitResult::SetParameterNames(const std::vector<std::string>& names_){
    fParameterNames = names_;
}

std::vector<std::string>
FitResult::GetParameterNames() const{
    return fParameterNames;
}

void
FitResult::Print() const{
    if(fParameterNames.size() != fBestFit.size())
        throw NotFoundError(Formatter() << "Expected one name for each parameter - got " 
                            << fParameterNames.size() << " names and " << fBestFit.size() << " params"
                            );

    std::cout << "Fit Result: " << std::endl;
    for(size_t i = 0; i < fParameterNames.size(); i++){
        std::cout << fParameterNames.at(i) << "\t" 
                  << fBestFit.at(i)
                  << std::endl;
    }        
}

void
FitResult::SaveAs(const std::string& fileName_) const{
  time_t t = time(0);
  struct tm * now = localtime(&t);
  
  std::ofstream fs(fileName_.c_str());
  fs << "-----------------------------------"
	 << "-----------------------------------" << std::endl
	 << "OXSX Fit Result: "
	 << now -> tm_year + 1900 << '-' 
	 << now -> tm_mon + 1 << '-'
	 << now -> tm_mday
	 << "  "
	 << now -> tm_hour
	 << ":"
	 << now -> tm_min
	 << ":"
	 << now -> tm_sec
	 << std::endl
	 << "-----------------------------------"
	 << "-----------------------------------" << std::endl
	 << std::endl;
	   
  fs << "Best Fit Values: " << std::endl << std::endl;
  for(size_t i = 0; i < fParameterNames.size(); i++){
	fs << std::setw(25) 
	   << fParameterNames.at(i) << "\t\t" 
	   << std::setw(10) 
	   << fBestFit.at(i)
	   << std::endl;
  }        	  
  fs.close();
}
