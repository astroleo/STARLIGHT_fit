from astropy.io import fits
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import numpy as np

##
## masking is problematic when re-binning to log space and especially when convolving with a large kernel to fit the SSP resolution (too many points may get masked); this routine will take all masked points and interpolate in between them linearly
## TODO: consider to increase the uncertainty of the "unmasked" points by some factor
##
## Wavelengths are in nm here!
def unmask_spectra(infile, qcplot=True, tac=False):
	hdu=fits.open(infile)
	t=hdu[1].data
	w=t['WAVE']
	f=t['FLUX']
	n=t['NOISE']
	q=t['QUAL']

	##
	## tac keyword specifies whether spectrum has been telluric absorption corrected
	##    using molecfit; in this case there is an additional quality flag to take into account
	if tac==True:	
		tac_q=t['tacqual']
		q=np.all([q,tac_q],axis=0)

	##
	## interpolation function
	interp_fct=interp1d(w[q==1],f[q==1])
	interp_fct_n=interp1d(w[q==1],n[q==1])
	##
	## truncate wave and qual array since at the edges there are often only masked data
	## which can not be correctly interpolated -- hardcoded for VIS range!
	croprange=[545,1015]
	## try without cropping
#	croprange=[0,1e7]

	cropmask=(w>croprange[0]) & (w<croprange[1])
	wcrop=w[cropmask]
	qcrop=q[cropmask]
	f_interpolated = interp_fct(wcrop[qcrop!=1])
	n_interpolated = interp_fct_n(wcrop[qcrop!=1])
	##
	## replace bad values by interpolated values
	f_clean=f.copy()
	n_clean=n.copy()
	qcropmask=(w>croprange[0]) & (w<croprange[1]) & (q!=1)
	f_clean[qcropmask]=f_interpolated
	n_clean[qcropmask]=n_interpolated
	##
	## QC plot
	if qcplot == True:
		plt.plot(w,f,'grey')
		plt.plot(w[q!=1],f[q!=1],'rx')
		plt.plot(wcrop[qcrop!=1],f_interpolated,'g.')
		plt.plot(w,f_clean,'b')

		plt.plot(w,n,'grey')
		plt.plot(w[q!=1],n[q!=1],'rx')
		plt.plot(wcrop[qcrop!=1],n_interpolated,'g.')
		plt.plot(w,n_clean,'purple')

		#plt.xlim([830,890])
		m=np.median(f)
		s=np.std(f)
		plt.ylim([0,m+s])
		outfile_qc=infile.split(".fits")[0]+'_unmask_QC.pdf'
		plt.savefig(outfile_qc)
		print("Saved QC plot to ", outfile_qc)
		plt.close()
	##
	## end of QC plot
	return(w,f_clean,n_clean)