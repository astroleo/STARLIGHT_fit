import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import os
from unmask_spectra import unmask_spectra
from matplotlib import pyplot as plt

import pdb

##
## arguments:
##    calspec_file   calibrated spectrum of science object
##    outfile        path of output (txt) file
##    z              redshift of object
##    wrange         restframe wavelength range to use (must obviously be smaller 
##                     than de-redshifted wavelength range of input spectrum)
##                     wrange is given in Angstrom, inclusive range, i.e. start *and* end points are included
##
def write_spec_starlight(calspec_file,outfile,z,wrange):
	wave,flux,noise = unmask_spectra(calspec_file)
	wave*=10 ## in Angstrom
	
	##
	## de-redshift spectrum
	wave/=(1+z)
	
	##
	## check that input wavelength range is larger than output wavelength range
	if (wave[0] > wrange[0]) or (wave[-1] < wrange[-1]):
		raise ValueError("input wavelength range must be larger than output wavelength range")
	##
	## normalize flux and noise
	m=np.median(flux)
	flux/=m
	noise/=m
	
	##
	## define output wavelength vector
	wsampling=1 ## wavelength sampling in Angstrom
	lamgrid=wrange[0] + wsampling * np.arange(1+np.floor((wrange[1]-wrange[0])/wsampling))

	##
	## convolve to resolution of spectral base, i.e. BC03 res = 3 \AA
	##
	## XSHOO resolution is 12600 (VIS) independent of lambda, i.e. choose kernel 
	##    appropriate for central wavelength
	xres=12600
	FWHM_gal=np.average(wrange)/xres
	FWHM_tem=3
	FWHM_dif = np.sqrt(FWHM_tem**2 - FWHM_gal**2)
	sigma = FWHM_dif/2.355/(wave[2]-wave[1]) # Sigma difference in pixels
	
	flux=gaussian_filter1d(flux,sigma)
	##
	## interpolation is flux-conserving when accounting for different spectral pixel size
	##    dl is 0.2 in UVB and VIS, 0.6 in NIR
	##
	f=interp1d(wave,flux)
	iflux=f(lamgrid)
	
	##
	## interpolate inverse variance (since this is the additive quantity)
	##
	## interpolation is not enough; need to also reduce the noise!
	v=interp1d(wave,1/(noise**2))
	ivariance=v(lamgrid)
	inoise=1/np.sqrt(ivariance)
	##
	## reduce noise by 1/sqrt(N) where N is number of bins averaged over
#	pdb.set_trace()
	N = wsampling/np.diff(wave)
	N2=np.hstack([N[0],N]) ## np.diff reduces size of array by 1
	Ni = interp1d(wave,N2)
	Ninterp = Ni(lamgrid)

	fac=np.sqrt(Ninterp)
	inoise/=fac

	## masks for transition regions of XSHOOTER spectral arms
	imask=np.zeros(iflux.shape)
#	add_mask_telluric(lamgrid, imask, [5530,5600], z)
#	add_mask_telluric(lamgrid, imask, [10100,10200], z)

	## Quality Check plot
	##
	plt.plot(wave,flux,label="intrinsic resolution")
	plt.plot(lamgrid,iflux,label="downsampled resolution")
	##
	## pick region around CaT if within region
	CaT1=8400
	CaT2=8800
	
	if (8000 > wrange[0]) & (9000 < wrange[1]):
		plt.xlim([8000,9000])
		plt.ylim([0.8,1.2])
	else:
		plt.xlim(wrange)
		plt.ylim([0,2])

	plt.legend(loc=2)
	outfile_qc=outfile.split('.txt')[0]+'_interpolate_QC.pdf'
	plt.savefig(outfile_qc)
	print("Saved QC plot to ", outfile_qc)
	
	np.savetxt(outfile,np.transpose((lamgrid,iflux,inoise,imask)),fmt=('%5.0f.   %8.3f   %8.3f    %3.0f'))
	print("Wrote ", outfile)