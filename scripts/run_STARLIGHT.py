from matplotlib import pyplot as plt
import numpy as np
import os
import subprocess
import glob
from SLplotlib import *

##
## helper function to parse STARLIGHT config files and adapt plotting range etc.
##
## currently only gets config file for first fit
##
def get_SL_params(infile):
	a=subprocess.run("grep Olsyn_ini " + infile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
	Olsyn_ini=np.float(a.stdout.strip())
	a=subprocess.run("grep Olsyn_fin " + infile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
	Olsyn_fin=np.float(a.stdout.strip())
	a=subprocess.run("head -n 16 " + infile + " | tail -n 1 | awk '{print $2}'", shell=True, universal_newlines=True, stdout=subprocess.PIPE)
	configfile=a.stdout.strip()

	a=subprocess.run("grep llow_norm " + configfile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
	llow_norm=np.float(a.stdout.strip())
	a=subprocess.run("grep lupp_norm " + configfile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
	lupp_norm=np.float(a.stdout.strip())
	return(Olsyn_ini,Olsyn_fin,llow_norm,lupp_norm)
 
dir=os.getenv('HOME')+'/STARLIGHT/'
os.chdir(dir)
infile="XSHOO.in"
subprocess.run("./StarlightChains_v04.exe < " + infile, shell=True)
Olsyn_ini,Olsyn_fin,llow_norm,lupp_norm = get_SL_params(infile)

dir=os.getenv('HOME')+'/STARLIGHT/out_dir/'
os.chdir(dir)
subprocess.call("../scripts/extract_results.sh",shell=True)

i=1
for file_syn in glob.glob("*.synspec"):
	if i > 1:
		continue
	id=file_syn.split(".")[0]
	plotSL(id,pdf=True, Olsyn_ini=Olsyn_ini, Olsyn_fin=Olsyn_fin, llow_norm=llow_norm, lupp_norm=lupp_norm)
	i+=1

#plot_popvec_hist()