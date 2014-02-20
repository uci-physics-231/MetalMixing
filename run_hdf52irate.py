import sys
import imp
import numpy as np
import Ngadget as gadget

# run run_hdf52irate.py <icnum> <runnum> <snap1> <snap2>
# run run_hdf52irate.py G00 0 1

workpath = "/Users/coralwheeler/Work/Dwarf/MetalMixing"

#######################
if __name__ == "__main__":
    if len(sys.argv)<3:
        print "Usage: python run_hdf52irate.py <icnum> <runnum> <snap>"
        sys.exit(1337)
    icnum=sys.argv[1] #initial conditions
    runnum=sys.argv[2] # gadget run
    snap1=int(sys.argv[3]) # snap number must not have any leading zeroes
    snap2=int(sys.argv[4])

s8="0.801";ns="0.963";omegaB="0.0449"
lunits=np.array([3.08568025e21,-1,1])

filein="%s/%s_%s_%03d.hdf5"%(workpath,icnum,runnum,int(snap1))
fileout="%s/%s_%s_%03d-irateTEST.hdf5"%(workpath,icnum,runnum,int(snap1))

print filein
print fileout

gadget.gadget_hdf5_to_irate(filein,fileout,snap2,lunits=lunits,s8=s8,ns=ns,omegaB=omegaB)
