from matplotlib import pyplot as plt
import numpy as np
import os
import glob
from SLplotlib import *

dir=os.getenv('HOME')+'/STARLIGHT/out_dir/'
os.chdir(dir)

i=1
for file_syn in glob.glob("*.synspec"):
#	if i > 1:
#		continue
	id=file_syn.split(".")[0]
	plotSL(id,pdf=True)
	i+=1

#plot_popvec_hist()