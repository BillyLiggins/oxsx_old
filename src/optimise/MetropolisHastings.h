#ifndef __OXSX_METROPOLIS_HASTINGS__
#define __OXSX_METROPOLIS_HASTINGS__
#include <Optimiser.h>
#include <FitResult.h>
#include <Histogram.h>

class TestStatistic;
class MetropolisHastings : public Optimiser{
 public:
     MetropolisHastings() : fBurnIn(3000), 
                            fMaxIter(100000), fThinFactor(1), 
                            fMaxVal(0), fFlipSign(false), 
                            fTestStatLogged(false), pTestStatistic(NULL)
                            {}               
    
    const FitResult& Optimise(TestStatistic*); 

    unsigned GetBurnIn() const;
    void     SetBurnIn(unsigned);
    
    unsigned GetThinFactor() const;
    void     SetThinFactor(unsigned);
    
    unsigned GetMaxIter() const;
    void     SetMaxIter(unsigned);

    double   GetRejectionRate() const;

    std::vector<double> GetMaxima() const;
    void SetMaxima(const std::vector<double>&);

    std::vector<double> GetMinima() const;
    void SetMinima(const std::vector<double>&);

    const std::vector<double>& GetSigmas() const;
    void   SetSigmas(const std::vector<double>&);

    bool GetFlipSign() const;
    void SetFlipSign(bool);
    
    bool GetTestStatLogged() const;
    void SetTestStatLogged(bool b_);

    void SetInitialTrial(const std::vector<double>&);
    std::vector<double> GetInitialTrial() const;

 private:
    unsigned fBurnIn;
    unsigned fThinFactor;
    unsigned fMaxIter;
    unsigned fNDims;
    bool     fTestStatLogged;
    
    TestStatistic* pTestStatistic;
    std::vector<double> fMaxima;
    std::vector<double> fMinima;
    std::vector<double> fSigmas;
    std::vector<double> fInitialTrial;

    bool fFlipSign;

    double fRejectionRate;

    std::vector< std::vector<double> > fSample;
    Histogram fHist;
    FitResult fFitResult;
    std::vector<double> fBestFit;
    double fMaxVal;
    
    std::vector<double> JumpDraw(const std::vector<double>& thisStep_) const; 
    inline double JumpProbRatio(const std::vector<double>& thisStep_, 
                                const std::vector<double>& proposedStep_) const {return 1;}

    bool   StepAccepted(const std::vector<double>& thisStep_,
                        const std::vector<double>& proposedStep_);
    


};
#endif

