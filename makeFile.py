#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
No doubt some of the worst code you have ever written.

If you run this then prune_mod.py is called to extract branches out of the files listed below.
"""
import os 
import glob 
import numpy as np
import pprint as pp


def massPrune():
    _Bi210PATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedBi210/TeLoadedBi210_r402_s0_p0.ntuple.root")
    _Po210PATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedPo210/TeLoadedPo210_r1_s0_p3.ntuple.root")
    _C14PATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedC14/TeLoadedC14_r10_s0_p2.ntuple.root")
    _2n2bPATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedTe130_0n2b_ntuples/TeLoadedTe130_0n2b_r100_s0_p1.ntuple.root")
    _0n2bPATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedTe130_2n2b/TeLoadedTe130_2n2b_r*")
    # Tl208LIST= ls /data/snoplus/OfficialProcessing/production_5_0/

    paths= [[_Bi210PATH, "Bi210"],[_Po210PATH,"Po210"],[_C14PATH, "C14"],[_2n2bPATH,"2n2b"],[_0n2bPATH,"0n2b"]]
    branches = ["energy","nhits"]
    branchstr= " ".join(branches)
    print paths

    for path, name in paths:
        print path, " ", name
        os.system("python prune_mod.py -filepath "+path+"/\* -folder /data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/"+str(name)+"/ -branches "+branchstr)
        print("python prune_mod.py -filepath "+path+"/ -folder /data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/"+str(name)+"/ -branches "+branchstr)

def massHadd():
    
    _Bi210PATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedBi210/TeLoadedBi210_r402_s0_p0.ntuple.root")
    _Po210PATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedPo210/TeLoadedPo210_r1_s0_p3.ntuple.root")
    _C14PATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedC14/TeLoadedC14_r10_s0_p2.ntuple.root")
    _2n2bPATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedTe130_0n2b_ntuples/TeLoadedTe130_0n2b_r100_s0_p1.ntuple.root")
    _0n2bPATH= os.path.dirname("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedTe130_2n2b/TeLoadedTe130_2n2b_r*")
    # Tl208LIST= ls /data/snoplus/OfficialProcessing/production_5_0/

    paths= [[_Bi210PATH, "Bi210"],[_Po210PATH,"Po210"],[_C14PATH, "C14"],[_2n2bPATH,"2n2b"],[_0n2bPATH,"0n2b"]]
    pathnp=np.array(paths)
    pathList=["/data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/"+str(i) for i in pathnp[:,1]]

    filelist= [ glob.glob(i+"/*") for i in pathList]

    for j, files in zip(pathnp[:,1],filelist):
        filestr=" ".join(files)
        pp.pprint("hadd /data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/"+str(j)+"/complete_"+str(j)+".ntuple_oxsx.root "+filestr )
        pp.pprint( "++++++++++++++++++++++++++++++++++")
        print "Make sure you have sourced root (rat)"
        os.system("hadd -f /data/snoplus/liggins/year1/fitting/oxsx/testData/TeLoaded/"+str(j)+"/complete_"+str(j)+".ntuple_oxsx.root "+filestr )

if __name__ == '__main__':
    """ As said above this is terrible code and should be rewritten completely,
    however to use this just add the path and name to the paths varible. The
    path is the path to the folder which contains all the SNO+ ntuples and
    "name" is the name of the target folder e.g. testData/TeLoaded/Bi210/ Bi210
    is the name."""

    # massPrune()

    massHadd()
