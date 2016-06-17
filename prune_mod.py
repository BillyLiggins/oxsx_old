#! /usr/bin/env python

################################################################
# Prune away branches of a flat root tree to produce an ntuple #
#
# example usage:
# python prune_mod.py -b energy time -folder testData/solar/ -filepath "/data/snoplus/OfficialProcessing/production_5_3_0/Bi210/SolarBi210_*"
################################################################
import argparse
import os
import sys
import glob
from array import array
try:
    from ROOT import TFile, TTree, TNtuple, TChain
except ImportError:
    print "Cant find root lib is it in $PYTHONPATH?"
    sys.exit()


def make_ntup(file_name, tree_name,  branches, outfile, n_events, new_tree_name):
    if new_tree_name == "":
        new_tree_name = tree_name

    print file_name
    
    # Get the event tree 
    tree = TChain(tree_name)
    tree.Add(file_name)
    if not tree:
        print "Error: No tree named %s in %s" %(tree_name, file_name)
        sys.exit()
    
    # Check branches exist
    branches_avail = [x.GetName() for x in tree.GetListOfBranches()]
    for b in branches:
        if not b in branches_avail:
            print "Error branch '%s' not a branch in input tree" %(b)
            print "Branches available are: \n"
            print "\t".join(branches_avail)
            sys.exit()

    # output
    out_file = TFile(outfile, "RECREATE")
    nt       = TNtuple(new_tree_name, "", ":".join(branches))

    if(n_events < 0):
        n_events = tree.GetEntries()

    # loop over events and fill the branches of new ntuple
    for index, entry in enumerate(tree):
        if index > n_events:
            break
        vals = array('f', [entry.__getattr__(b) for b in branches])
        nt.Fill(vals)
        
        if (index % 100000 == 0):
            print index, "/", n_events
    # Save
    out_file.cd()
    nt.Write()
    out_file.Close()

    print "Written %i entries of branch(es) '%s' \nto tree %s  \nin file %s" %(n_events, 
                                                                             ":".join(branches), 
                                                                             new_tree_name, outfile)
def convertFiles(filelist, tree_name,  branches, outfile, n_events, new_tree_name,folder):
    for i in filelist:
        filename__=os.path.split(i)[1]

        outfile = os.path.splitext(filename__)[0]
        outfile += "_oxsx.root"
        outfile = folder+"/"+outfile
        outfile = outfile
        print "this is the outfile : ",outfile
        make_ntup(i, tree_name, branches,outfile, n_events, new_tree_name)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-filepath', type=str, default = "")
    parser.add_argument('-treename', metavar='-t', type=str, default = "output")
    parser.add_argument('-newtreename', type=str, default = "")
    parser.add_argument('-branches',metavar='-b', nargs="+", type=str)
    parser.add_argument('-outfile', metavar='-o', type=str, default = "")
    parser.add_argument('-nevents', metavar='-nev', type=int, default = -1)
    parser.add_argument('-folder', metavar='-f', type=str, default = "")
    args = parser.parse_args()
    
    # filelist=glob.glob("/data/snoplus/OfficialProcessing/production_5_0/TeLoadedTe130_2n2b/TeLoadedTe130_2n2b_*")
    # filelist=glob.glob("/home/billy/workspace/PhD/testData/Po210/*")
    filelist=glob.glob(args.filepath)
    print "+++++++++++++++++++FileList+++++++++++++++++++++++"
    print filelist
    print "you have to give this a filelist, branches and a folder!"

    convertFiles(filelist, args.treename, args.branches,args.outfile, args.nevents, args.newtreename,args.folder)
    

        
