from astroquery.simbad import Simbad
from scipy import constants
from write_spec_starlight import write_spec_starlight
import os

def run_starlight(object_id,object_name,ob_id):
	SLdir=os.getenv("HOME")+"/STARLIGHT"
	scifile=SLdir+"/spectra/"+object_id+"_"+str(ob_id)+".fits"
	outfile=SLdir+"/spectra/"+object_id+"_"+str(ob_id)+".txt"
	##
	## get redshift of object
	Simbad.add_votable_fields("rv_value")
	q=Simbad.query_object(object_name)
	rv=q['RV_VALUE'].quantity.data[0]*1000 ## rv in m/s
	z=rv/constants.c
	print("Redshift of {object_id} is {z:7.5f}".format(object_id=object_id, z=z))
	##
	## get spectrum from FITS file, interpolate bad pixels, convolve to BC03 resolution for STARLIGHT
	write_spec_starlight(scifile,outfile,z,[6000,9500])

object_id="ESO208"
object_name="ESO208-G021"
ob_id=1
run_starlight(object_id,object_name,ob_id)