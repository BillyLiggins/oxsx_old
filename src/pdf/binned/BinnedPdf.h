/**************************************************************************************************/
/* A BinnedPdf is an ND histogram along with a data representation that tells the pdf where it    */
/* looks for the data to pass to the histogram inside. Most methods just defer to the histogram.  */
/* The two are serparate classes, because other things are histograms e.g. LH space but not pdfs  */
/* in the sense that they have a data representation and examine events                           */
/**************************************************************************************************/

#ifndef __OXSX_BINNED_PDF__
#define __OXSX_BINNED_PDF__
#include <Histogram.h>
#include <Pdf.h>
#include <vector>

class BinnedPdf : public Pdf{
 public:
    BinnedPdf() {}
    BinnedPdf(const AxisCollection& axes_);
    BinnedPdf(const Histogram& histo_);
    Pdf*   Clone() const; 

    double operator() (const std::vector<double>& vals_) const;
    double Probability(const EventData&) const;
    double Integral()  const;
    void   Normalise();
    void   Scale(double s_);
    
    void   Fill(const std::vector<double>& vals_, double weight_ = 1);
    void   Fill(const EventData& data_, double weight_ = 1);
    void   Fill(double val_, double weight_ = 1);

    size_t FindBin(const std::vector<double>& vals_) const;
    size_t FindBin(const EventData& data_) const;

    std::vector<size_t> UnpackIndices(size_t bin_) const;
    size_t FlattenIndices(const std::vector<size_t>& indices_) const;

    const AxisCollection& GetAxes() const;
    void  SetAxes(const AxisCollection& axes_);

    std::vector<double> GetBinContents() const; 
    void SetBinContents(const std::vector<double>& data_);     

    std::vector<double> Means() const;
    std::vector<double> Variances() const;

    const Histogram& GetHistogram() const;
    void  SetHistogram(const Histogram&);

    double GetBinContent(size_t bin_) const;
    size_t GetNBins() const;

    void     AddBinContent(size_t bin_, double content_);
    void     SetBinContent(size_t bin_, double content_);
    void     Empty();
    unsigned GetNDims() const;
        
    BinnedPdf Marginalise(const std::vector<size_t>& indices_) const;
    BinnedPdf Marginalise(size_t index_) const;

    void SetDataRep(const DataRepresentation&);
    DataRepresentation GetDataRep() const;

 protected:
    DataRepresentation fDataRep;
    Histogram fHistogram;
};
#endif
