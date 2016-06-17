#include <ComponentManager.h>
#include <Exceptions.h>
#include <algorithm>
#include <iostream>

void 
ComponentManager::AddComponent(FitComponent*  componentPtr_){
    componentPtr_ -> MakeFittable();
    fComponents.push_back(componentPtr_);
    fParamCounts.push_back(componentPtr_ -> GetParameterCount());
    fTotalParamCount += componentPtr_    -> GetParameterCount();
    fComponentCount++;
}

void 
ComponentManager::SetParameters(const std::vector<double>& params_){
    if(params_.size() != fTotalParamCount)
        throw ParameterCountError("Component Manager", fTotalParamCount, 
                                  params_.size());

    if(!params_.size())
        return;

    std::vector<double>::const_iterator it = params_.begin();
    fComponents[0] -> SetParameters(std::vector<double>(it, 
                                                        it + fParamCounts.at(0)));
    for(size_t i = 1; i < fComponents.size(); i++){
        fComponents[i] -> SetParameters(std::vector<double>(it + fParamCounts.at(i-1),
                                                            it  + fParamCounts.at(i-1) +
                                                            + fParamCounts.at(i)
                                                            )
                                        );
    }
}

std::vector<double> 
ComponentManager::GetParameters() const{
    std::vector<double> params;
    for(size_t i = 0; i < fComponents.size(); i++){
        const std::vector<double>& comps = fComponents.at(i) -> GetParameters();
        params.insert(params.end(), comps.begin(), comps.end());
    }
    return params;
}

std::vector<std::string>
ComponentManager::GetParameterNames() const{
    std::vector<std::string> paramNames;
    for(size_t i = 0; i < fComponents.size(); i++){
        const std::vector<std::string>& comps = fComponents.at(i) -> GetParameterNames();
        paramNames.insert(paramNames.end(), comps.begin(), comps.end());
    }
    return paramNames;
}

int
ComponentManager::GetTotalParameterCount() const{
    return fTotalParamCount;
}

void
ComponentManager::Clear(){
    fTotalParamCount = 0;
    fComponentCount  = 0;
    fComponents.clear();
    fParamCounts.clear();
}

double
ComponentManager::GetParameter(const std::string& paramName_) const{
    for(size_t i = 0; i < fComponents.size(); i++){
        FitComponent* component = fComponents[i];
        std::vector<std::string> names = component->GetParameterNames();

        std::vector<std::string>::iterator it = std::find(names.begin(), names.end(), 
                                                          paramName_);
        if(it != names.end())
            return component->GetParameters()[it - names.begin()];
    }
    throw NotFoundError(Formatter() << "ComponentManager::GetParameter "
                        << " requested non-existent parameter "
                        << paramName_
                        );
}
