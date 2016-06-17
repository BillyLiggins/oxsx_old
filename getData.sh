#!/bin/bash

echo making folders
mkdir testData/TeLoaded/
mkdir testData/TeLoaded/Bi210/
mkdir testData/TeLoaded/Po210/
mkdir testData/TeLoaded/C14/ 
mkdir testData/TeLoaded/0n2b/ 
mkdir testData/TeLoaded/2n2b/
mkdir testData/TeLoaded/Tl208/



echo copying data over

scp -r liggins@hep.ph.qmul.ac.uk:/data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/Bi210/complete_Bi210.ntuple_oxsx.root testData/TeLoaded/Bi210/
scp -r liggins@hep.ph.qmul.ac.uk:/data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/Po210/complete_Po210.ntuple_oxsx.root testData/TeLoaded/Po210/
scp -r liggins@hep.ph.qmul.ac.uk:/data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/C14/complete_C14.ntuple_oxsx.root testData/TeLoaded/C14/
scp -r liggins@hep.ph.qmul.ac.uk:/data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/0n2b/complete_0n2b.ntuple_oxsx.root testData/TeLoaded/0n2b/
scp -r liggins@hep.ph.qmul.ac.uk:/data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/2n2b/complete_2n2b.ntuple_oxsx.root testData/TeLoaded/2n2b/


echo data copied
