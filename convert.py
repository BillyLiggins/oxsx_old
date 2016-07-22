import hadd_mod as hd
import os


extensions= ["Bi210","Bi212","Bi214","Po210","Po212"]


for i in extensions:
    if os.path.isfile("testData/Bab_temp/"+str(i)+"/Bab_"+str(i)+"_complete.ntuple_oxsx.root"):
        os.remove("testData/Bab_temp/"+str(i)+"/Bab_"+str(i)+"_complete.ntuple_oxsx.root")
    # print "testData/Bab_temp/"+str(i)+"/","testData/Bab_temp/"+str(i),"Bab_"+str(i)+"_complete.ntuple_oxsx.root",4
    hd.HADD("testData/Bab_temp/"+str(i)+"/","testData/Bab_temp/"+str(i),"Bab_"+str(i)+"_complete.ntuple_oxsx.root",4)

