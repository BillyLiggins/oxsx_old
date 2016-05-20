#include <Ensemble.h>
#include <TestStatistic.h>
#include <Exceptions.h>
#include <Optimiser.h>
#include <FitResult.h>

void
Ensemble::RunPseudoExperiment(TestStatistic* stat_){
    if(!fNParams){
        fNParams = stat_->GetParameterCount();

    }
    
    if(fNParams != stat_->GetParameterCount())
        throw ParameterCountError("Ensemble", fNParams, stat_->GetParameterCount(), " the expected number was set the first time you called RunPseudoExperiment");

    const FitResult& fitResult = fOptimiser -> Optimise(stat_);
    const std::vector<double>& pointEstimate = fitResult.GetBestFit();

}
