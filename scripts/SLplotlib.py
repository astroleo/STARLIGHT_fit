from matplotlib import pyplot as plt
from astropy.io import ascii
import numpy as np
import os
import glob

import pdb


##
## wishlist
##    - widget button to show plots of next object
##    - read llow_norm, lupp_norm from StVc04.C11.config, instead of hardcoding it here
##

##
## synRange: restrict plot to wavelength range in which the synthesis was done

def plotSL(id,pdf=False,synRange=True,plotSFH=True,llow_norm=6810,lupp_norm=6870,Olsyn_ini=3800,Olsyn_fin=10000):
	os.chdir(os.getenv("HOME")+"/STARLIGHT")
	
	file_obs="spectra/"+id+".txt"
	file_syn="out_dir/"+id+".out.synspec"
	file_masks="out_dir/"+id+".out.masks"
	file_popvec="out_dir/"+id+".out.popvec"

	spec_obs=ascii.read(file_obs, names=('wave','flux','noise','flag'))
	spec_syn=ascii.read(file_syn, names=('wave','fobs','fsyn','weight'))
	masks=ascii.read(file_masks, names=('lam1','lam2'))

	##
	## normalize spectrum the same way STARLIGHT does
	norm=np.median(spec_obs['flux'][(spec_obs['wave'] >= llow_norm) & (spec_obs['wave'] <= lupp_norm)])
	wobs=spec_obs['wave']
	fobs=spec_obs['flux']/norm
	nobs=spec_obs['noise']/norm
	flagobs=spec_obs['flag']

	wsyn=spec_syn['wave']
	fsyn=spec_syn['fsyn']


	###
	## wavelength range in which the fit was done
	###
	if synRange == True:
		Xcfg="XSHOO.in"
		w_min=0
		w_max=0
		for line in open(Xcfg):
			if (w_min != 0) & (w_max != 0):
				break
			a=line.split()
			if a[1] == "[Olsyn_ini]":
				w_min=np.float(a[0])
			if a[1] == "[Olsyn_fin]":
				w_max=np.float(a[0])
	else:
		w_min=wobs[0]
		w_max=wobs[-1]

#	print("w_min and w_max are: {w_min},{w_max}".format(w_min=w_min,w_max=w_max))

	##
	## later: read also Mask file to plot masked emission lines
	#with open('Mask.ESO093.BN','r') as f:
	#	for line in f:
	#		if i == 1:
	#			nlines=line
	#		i+=1

	if plotSFH == False:
		gridsize=(4,3)
		gridpos_spec=(0,0)
		gridpos_res=(3,0)
	if plotSFH == True:
		gridsize=(4,4)
		gridpos_spec=(0,0)
		gridpos_res=(3,0)
		gridpos_sfh=(2,3)

	plt.subplot2grid(gridsize,gridpos_spec,rowspan=3,colspan=3)
	plt.minorticks_on()
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.title(id)

	lw=0.5
	plt.plot(wobs, fobs, 'k-', linewidth=lw)
	
	##
	## plot observed spectrum in yellow where flagged
#	plt.plot(wobs[np.where((flagobs != 0) & (wobs < 7200))], fobs[np.where((flagobs != 0) & (wobs < 7200))], 'y|')
#	plt.plot(wobs[np.where((flagobs != 0) & (wobs > 7200))], fobs[np.where((flagobs != 0) & (wobs > 7200))], 'y|')
#	plt.plot(wobs[np.where(flagobs != 0)], fobs[np.where(flagobs != 0)], 'y|', linewidth=lw)
	
	ylim_lo = 0.5
	ylim_hi = 3
	
	w_min=Olsyn_ini
	w_max=Olsyn_fin
	
	plt.xlim((w_min,w_max))
#	plt.ylim((ylim_lo,ylim_hi))
	plt.ylim((0,5))
	plt.ylabel(r'$F_{\lambda}$' + ' [normalized]',fontsize=10)

	plt.plot(wobs, nobs, 'g', linewidth=lw)
	plt.plot(wsyn, fsyn, 'b', linewidth=lw)

	##
	## residual plot
	try:
		assert(fsyn.size < fobs.size)
	except:
		print("fsyn.size is greater than fobs.size")

	## bring fsyn and fobs to same length
#	pdb.set_trace()
	fres = fobs[(wobs >= wsyn[0]) & (wobs <= wsyn[-1])] - fsyn

	plt.subplot2grid(gridsize,gridpos_res,rowspan=1,colspan=3)
	plt.minorticks_on()
	plt.tick_params(axis='both', which='major', labelsize=8)

	plt.plot(wsyn, fres, 'k', linewidth=lw)
	##
	## mark points that have not been used for the fit
	##
	ix=np.where(spec_syn['weight'] <= 0)
	ms=1.0
	plt.plot(wsyn[ix], fres[ix], 'rx', markersize=ms)

	plt.xlim((w_min,w_max))
	plt.ylim((-0.3,0.5))
	plt.xlabel(r'$\lambda$' + ' [$\AA$]',fontsize=10)
	plt.ylabel(r'Residual',fontsize=10)
	
	for m in masks:
		plt.plot(wsyn[(wsyn >= m['lam1']) & (wsyn <= m['lam2'])], 
			fres[(wsyn >= m['lam1']) & (wsyn <= m['lam2'])], 'r', linewidth=lw)
	
	if plotSFH == True:
		plt.subplot2grid(gridsize,gridpos_sfh,rowspan=2,colspan=1)
		plot_sfh(file_popvec)
		plt.subplots_adjust(left=0.08,bottom=0.09,right=0.96,top=0.93,wspace=0.55,hspace=0.32)


	if pdf == True:
		file_plot="plot/"+id+".pdf"
		plt.savefig(file_plot)
		plt.tight_layout()
		plt.close()
	else:
		plt.show()

###
## get population vector as functions of light and mass
##
## nice to have:
##    - more generic (e.g. automatically determine no. of metallicities from output data)
##    - check whether age distributions are really the same for all metallicities
##    - distinguish light/mass weighted
###
def plot_sfh(file_popvec):
	popvec=np.genfromtxt(file_popvec)

	## generate popvecs (light/mass weighted) for the three metallicities
	popvec_sub=popvec[np.where(popvec[:,5]==0.004)]
	popvec_sol=popvec[np.where(popvec[:,5]==0.02)]
	popvec_sup=popvec[np.where(popvec[:,5]==0.05)]

	popvec_lw_sub=popvec_sub[:,1]
	popvec_lw_sol=popvec_sol[:,1]
	popvec_lw_sup=popvec_sup[:,1]

	agevec_log=np.log10(popvec_sub[:,4])

	plt.bar(agevec_log,popvec_lw_sub, width=0.15,color='b',label='sub-solar')
	plt.bar(agevec_log,popvec_lw_sol, width=0.15,color='g',bottom=popvec_lw_sub,label='solar')
	plt.bar(agevec_log,popvec_lw_sup, width=0.15,color='r',bottom=popvec_lw_sol+popvec_lw_sub,label='super-solar')
	
	plt.xticks((6,7,8,9,10),('6','7','8','9','10'))


	plt.legend(fontsize='xx-small')

	plt.xlabel("Log age [yr]")
	plt.ylabel("x_j (light) [%]")


def get_sfh_summary(file_popvec):
	popvec=np.genfromtxt(file_popvec)
	y=np.where(popvec[:,4] <= 2.5e7)
	o=np.where(popvec[:,4] > 1.4e9)
	i=np.where((popvec[:,4] > 2.5e7) & (popvec[:,4] <= 1.4e9))
	x_y=np.sum(popvec[y,1])
	x_i=np.sum(popvec[i,1])
	x_o=np.sum(popvec[o,1])
	##
	## normalize
	x_total=0.01*(x_y+x_i+x_o)
	x_y/=x_total
	x_i/=x_total
	x_o/=x_total
	
	return x_y, x_i, x_o

def plot_popvec_hist():
	dir=os.getenv('HOME')+'/STARLIGHT/out_dir/'
	os.chdir(dir)
	
	x_y_active=[]
	x_i_active=[]
	x_o_active=[]
	x_y_control=[]
	x_i_control=[]
	x_o_control=[]

	for file_syn in glob.glob("*.synspec"):
		id=file_syn.split(".")[0]
	
		file_popvec=id+".out.popvec"
		x_y, x_i, x_o = get_sfh_summary(file_popvec)
	
		if (id == "NGC4593") | (id == "MCG514") | (id == "NGC1365") | (id == "NGC2110") | (id == "NGC2992") | (id == "NGC3081"):
			print(id + " (AGN)")
			print(round(x_y), round(x_i), round(x_o))
			x_y_active.append(x_y)
			x_i_active.append(x_i)
			x_o_active.append(x_o)
		else:
			print(id + " (control galaxy)")
			print(round(x_y), round(x_i), round(x_o))
			x_y_control.append(x_y)
			x_i_control.append(x_i)
			x_o_control.append(x_o)

	plt.subplot(311)
	plt.hist((x_y_active,x_y_control),color=('b','r'),label=('AGN','control'))
	plt.legend()
	plt.xlim(0,100)
	plt.ylim(0,7)
	plt.title('Young population')
	plt.ylabel('Number of galaxies')

	plt.subplot(312)
	plt.hist((x_i_active,x_i_control),color=('b','r'),label=('AGN','control'))
	plt.xlim(0,100)
	plt.ylim(0,7)
	plt.title('Intermediate-age population')
	plt.ylabel('Number of galaxies')

	plt.subplot(313)
	plt.hist((x_o_active,x_o_control),color=('b','r'),label=('AGN','control'))
	plt.xlim(0,100)
	plt.ylim(0,7)
	plt.title('Old population')
	plt.xlabel('Light fraction')
	plt.ylabel('Number of galaxies')

	plt.tight_layout()

	plt.savefig('../plot/SFH_compare.pdf')
	plt.close()

if __name__ == "__main__":
    import sys
    plotSL(sys.argv[1])
