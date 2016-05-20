#ifndef __OXSX_ENSEMBLE__
#define __OXSX_ENSEMBLE__
#include <Histogram.h>

class TestStatistic;
class Optimiser;
class Ensemble{
 public:
    Ensemble(Optimiser* optimiser_) : fNParams(0), fOptimiser(optimiser_) {}
    void RunPseudoExperiment(TestStatistic*);

 private:
    Optimiser* fOptimiser;
    std::vector<Histogram> fFits;
    std::vector<Histogram> fBiases;
    std::vector<Histogram> fPulls;

    int fNParams;
};
#endif
