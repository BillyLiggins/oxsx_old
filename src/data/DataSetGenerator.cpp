#include <DataSetGenerator.h>
#include <DataSet.h>
#include <OXSXDataSet.h>
#include <Exceptions.h>
#include <math.h>
#include <EventData.h>
#include <stdlib.h>
#include <Rand.h>

void
DataSetGenerator::SetDataSets(const std::vector<DataSet*> sets_){
    fDataSets = sets_;
}

void
DataSetGenerator::SetExpectedRates(const std::vector<double>& rates_){
    fExpectedRates = rates_;
}

OXSXDataSet
DataSetGenerator::ExpectedRatesDataSet() const{
    if(fExpectedRates.size() != fDataSets.size())
        throw LogicError("Can't generate fake data: need one rate exactly for each data set");

    if(!fDataSets.size())
        throw LogicError("Can't generate fake data: no input data sets");

    OXSXDataSet dataSet;    
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
        unsigned expectedCounts = round(fExpectedRates.at(i));

        for(unsigned j = 0; j < expectedCounts;){
            EventData event_ = RandomEvent(i);
            
            // Add event to new data set if it passes cut
            if (fCuts.PassesCuts(event_)){
                dataSet.AddEntry(event_);
                j++;
            }
        }
    }
    
    return dataSet;
}

OXSXDataSet
DataSetGenerator::PoissonFluctuatedDataSet() const{
    if(fExpectedRates.size() != fDataSets.size())
        throw LogicError("Can't generate fake data: need one rate exactly for each data set");
    if(!fDataSets.size())
      throw LogicError("Can't generate fake data: no input data sets");

    OXSXDataSet dataSet;    
    dataSet.SetObservableNames(fDataSets.at(0)->GetObservableNames());
    for(size_t i = 0; i < fDataSets.size(); i++){
        int counts = Rand::Poisson(fExpectedRates.at(i));

        for(unsigned j = 0; j < counts;){
            EventData event_ = RandomEvent(i);
            
            // check it passes the cuts                            
            if (fCuts.PassesCuts(event_)){
                dataSet.AddEntry(event_);
                j++;
            }
        }
    }
        
    return dataSet;
}
    


EventData
DataSetGenerator::RandomEvent(size_t handleIndex_) const{
    unsigned eventNum = Rand::Shoot(fDataSets.at(handleIndex_)->GetNEntries());
    return fDataSets.at(handleIndex_)->GetEntry(eventNum);
}

void
DataSetGenerator::AddDataSet(DataSet* data_, double rate_){
    fDataSets.push_back(data_);
    fExpectedRates.push_back(rate_);
}

void
DataSetGenerator::SetCuts(const CutCollection& cuts_){
    fCuts = cuts_;
}

