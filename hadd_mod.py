#! /usr/bin/env python

################################################################
# This script is a simple hadd wrapper there are better online. 
#
# example usage:
# python hadd_mod.py -i "../testData/alpha/output/ntuple/" -n 4 -f testy.ntuple.root -o testData
#
################################################################

import argparse
import os
import sys
import glob
import numpy as np 
from array import array
try:
    from ROOT import TFile, TTree, TNtuple, TChain
except ImportError:
    print "Cant find root lib is it in $PYTHONPATH?"
    sys.exit()

def HADD(inputFolder,outputFolder,filename, numOfFiles):
    print inputFolder+"*"
    
    inputList=glob.glob(inputFolder+"*");
    np.random.shuffle(inputList)

    if (len(inputList)< numOfFiles):
        numOfFiles= len(inputList)

    # print inputList 
    # print numOfFiles

    print "hadd "+outputFolder+"/"+filename+" "+" ".join(inputList[:numOfFiles])
    os.system("hadd "+outputFolder+"/"+filename+" "+" ".join(inputList[:numOfFiles]))

if __name__ == "__main__":

    print os.getcwd()
    parser = argparse.ArgumentParser()
    parser.add_argument('-inputFolder',metavar='-i', type=str, default = "")
    parser.add_argument('-outputFolder',metavar='-o', type=str, default = ".")
    parser.add_argument('-filename',metavar='-f', type=str, default = "hadd_mod_output.nutple.root")
    parser.add_argument('-numOfFiles',metavar='-n', type=int, default = "")
    args = parser.parse_args()
    

    HADD(args.inputFolder,args.outputFolder,args.filename,args.numOfFiles)
    

        
