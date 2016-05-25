#ifndef __OXSX_FAKE_PDF_GENERATOR__
#define __OXSX_FAKE_PDF_GENERATOR__
#include <BinnedPdf.h>

class FakePdfGenerator{
 public:
    void SetPdfs(const std::vector<BinnedPdf>& pdfs_) {fPdfs = pdfs_;}
    std::vector<BinnedPdf> GetPdfs(const std::vector<BinnedPdf>& pdfs_) const {return fPdfs;}

    void SetRates(const std::vector<double> rates_) {fRates = rates_;}
    std::vector<double> GetRates() const {return fRates;}
    
    BinnedPdf ExpectedRatesPdf() const;
    BinnedPdf PoissonFluctuatedPdf() const;

 private:
    std::vector<BinnedPdf> fPdfs;
    std::vector<double>    fRates;
    size_t RandomBin(size_t pdfIndex_) const;
};
#endif
