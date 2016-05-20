#include <ParameterManager.h>
#include <DoubleParameter.h>
#include <Exceptions.h>
#include <iostream>

ParameterManager::~ParameterManager(){
    for(size_t i = 0; i < fParamPtrs.size(); i++)
        delete fParamPtrs[i];
}

void
ParameterManager::SetParameters(const std::vector<double>& vals_){
    if(vals_.size() != GetParameterCount())
        throw ParameterCountError("ParameterManager", 
                                  GetParameterCount(),
                                  vals_.size());

    for(size_t i = 0; i < vals_.size(); i++)
        fParamPtrs.at(i) -> Set(vals_.at(i));
}

std::vector<double>
ParameterManager::GetParameters() const{
    std::vector<double> params;
    params.reserve(fParamPtrs.size());

    for(size_t i = 0; i < fParamPtrs.size(); i++)
        params.push_back(fParamPtrs.at(i)->Get());
    return params;
}

std::vector<std::string>
ParameterManager::GetParameterNames() const{
    return fNames;
}

size_t
ParameterManager::GetParameterCount() const{
    return fParamPtrs.size();
}

void 
ParameterManager::Add(FitParameter* p_, const std::string& n_){
    fParamPtrs.push_back(p_);
    fNames.push_back(n_);
}

void 
ParameterManager::AddDouble(double& d_, const std::string& n_){
    fParamPtrs.push_back(new DoubleParameter(d_));
    fNames.push_back(n_);
}

void
ParameterManager::Clear(){
    fParamPtrs.clear();
    fNames.clear();
}
