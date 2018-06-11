# (Not-so-)Simple script to find GBNCC beams on zuul and, optionally, process with PRESTO     ###
# Processing requires environment set-up before script execution
# Can make plot of beams (observed, detected, and unobserved) as beamcheck_PSR.png
# Write out information about pulsars found as GBNCC_ATNF_beamcheck_[comp]_[n].txt
# Combine results from multiple computer systems with combine_beamcheck.py
# Make individual summaries with GBNCC_individual_summaries.py
# Make nice plots of profiles with GBNCC_profile_plots.py
# Written by Renee Spiewak
# Last edit Jun. 5, 2017

# -h,--help   Print help page with flags listed
# -pos        ra and dec in hh:mm:ss.ss (-)dd:mm:ss.ss form
# -p          period in s
# -pd         period derivative in s/s (default: 0)
# -dm         DM in pc/cc
# -f          file name (assumes file has columns of NUMBER, JNAME, RA (deg), DEC (deg), P0 (s), P1 (s/s), DM (cm^-3))
# -beam       beam number as integer
# -n          number of profile bins (default: 50)
# -nsub       number of channels (default: 128) 
# -npart      number of subintegrations (default: 40)
# -angle      max angular offset between beam and pulsar in deg (default: 0.5 deg)
# -snr        minimum S/N for detection (default: 6) 
# -proc       process files (default: no processing)
# -cen        process central (original) beam
# -comp       select computer system (only 'zuul' or 'GB' allowed)
# -pub        Only process published pulsars
# -rfi        Use Scott's options to remove RFI
# -nonew      Do not copy new fits files into directories 
# -out	      Choose output location 

#########################################################

import numpy as np
import sys, os
from glob import glob
import subprocess as sproc
#import matplotlib.pyplot as plt
from time import strftime
import argparse as ap
import astropy.coordinates as coord
from astropy.table import Table,Column
from astropy import units as u, constants as c
from astropy.coordinates import SkyCoord
# Function to find distance from Joe Swiggum

def ang_offset(lon1,lat1,lon2,lat2):
    d2r = np.pi/180.
    term1 = np.sin(lat1*d2r)*np.sin(lat2*d2r)
    term2 = np.cos(lat1*d2r)*np.cos(lat2*d2r)*np.cos(lon1*d2r-lon2*d2r)
    offset = np.arccos(term1+term2)
    return offset/d2r


def proc_args():
    """
    Intelligently read in arguments
    Check for errors

    Returns
    ---------
    dictionary containing relevant values

    """
    pars = ap.ArgumentParser(description="Locate and process data \
on known pulsars")
    pars.add_argument("--comp",action="store",choices=["zuul","GB"],
                      required=True,help="select computer system")
    pars.add_argument("--pos",action="store",
                      help="ra and dec ('hh:mm:ss.ss (-)dd:mm:ss.ss')")
    pars.add_argument("-p","--period",action="store",default=1.0,type=float,
                      help="period in s")
    pars.add_argument("-pd","--pdot",action="store",default=0,type=float,
                      help="period derivative in s/s")
    pars.add_argument("--dm",action="store",default=50,type=float,
                      help="DM in pc/cc")
    pars.add_argument("-f","--file",action="store",
                      help="file name (assumes file has columns of \
NUMBER, JNAME, RA, DEC, P0, P1, DM)")
    pars.add_argument("--beam",action="store",
                      help="beam number as integer")
    pars.add_argument("-n","--nbin",action="store",default=50,type=int,
                      help="number of profile bins")
    pars.add_argument("--nsub",action="store",default=128,type=int,
                      help="number of channels")
    pars.add_argument("--npart",action="store",default=40,type=int,
                      help="number of subintegrations")
    pars.add_argument("--angle",action="store",default=0.5,type=float,
                      help="max offset between beam and pulsar in deg")
    pars.add_argument("--snr",action="store",default=6,type=float,
                      help="minimum S/N for detection")
    pars.add_argument("--proc",action="store_true",
                      help="process files")
    pars.add_argument("--center",action="store_true",dest="find_psr",
                      help="process central (original) beam")
    pars.add_argument("--pub",action="store_true",
                      help="only process published pulsars")
    pars.add_argument("--rfi",action="store_true",
                      help="use Scott's options to remove RFI")
    pars.add_argument("--nonew",action="store_false",dest="pull_new",
                      help="do not copy new fits files into directories")
    pars.add_argument("-o","--out",action="store",help="output location (default is /users/rspiewak/pulsars/)")

    args = vars(pars.parse_args())

    if args['pos'] != None and len(args["pos"].split()) > 1:
        ra_astro = coord.Longitude(args["pos"].split()[0],unit="hourangle")
        dec_astro = coord.Latitude(args["pos"].split()[1],unit="deg")
        args["ra_psr"] = ra_astro.value
        args["dec_psr"] = dec_astro.value

    return args

# Begin full program 
args = proc_args()

p0_psr = args["period"]
pd_psr = args["pdot"]
dm_psr = args["dm"]
file_psr = args["file"]
beam_psr = args["beam"]
nbin = args["nbin"]
nsub = args["nsub"]
npart = args["npart"]
ang_max = args["angle"]
snr_min = args["snr"]
proc_psr = args["proc"]
find_psr = args["find_psr"]
use_comp = args["comp"]
f_pub = args["pub"]
rfi_fil = args["rfi"]
pull_new = args["pull_new"]
if "ra_psr" in args:
    ra_psr = args["ra_psr"]
    dec_psr = args["dec_psr"]


home_dir = "/users/rspiewak/pulsars/"
if "out" in args:
    output_loc = args["out"]
else:
    output_loc = home_dir
if use_comp == "zuul":
    work_dir = os.getcwd()+'/' #"/lustre/cv/projects/GBNCC/renee/"
    data_dir = "/lustre/cv/projects/GBNCC/"
elif use_comp == "GB":
    work_dir = "/home/scratch/rspiewak/beamcheck_all/"
    data_dir = "/lustre/pulsar/survey/AGBT09C_057/"

# Get info on all GBNCC past and future beams from Scott Ransom's file

beam_num = []
ra_astro = []
dec_astro = []

if not os.path.isfile('./GBNCC_pointings.fits'):
	filename='/users/sransom/GBNCC/GBNCC_posns_by_dec_ALLGBTSKY.txt'
	np.fromfile(filename)
	t=Table.read(filename,format='ascii',names=('pointing','RA','Dec'))
	s=SkyCoord(t['RA'],t['Dec'],unit=('hour','deg'))
	t.add_column(Column(s.ra,name='RAdeg'))
	t.add_column(Column(s.dec,name='Decdeg'))
	t.write('GBNCC_pointings.fits')
t=Table.read('GBNCC_pointings.fits')
ra_astro=t['RAdeg']
dec_astro=t['Decdeg']


beam_num = np.array(beam_num)
ra_full = np.array(ra_astro)
dec_full = np.array(dec_astro)
num_north=0

# Check file information (assumes order of data in file)
if not file_psr is None:
    num_all_psr = []
    name_all_psr = []
    ra_all_psr = []
    dec_all_psr = []
    p0_all_psr = []
    pd_all_psr = []
    dm_all_psr = []
    stat_all_psr = []
    read_all_psr = open(file_psr, 'r')
    pulsar_list = []
    for line in read_all_psr:
        if line.split()[1][0] == "J" or \
           (line.split()[1][0] == "C" and line.split()[2][0] != "P"):
            num_all_psr.append(int(line.split()[0]))
            name_all_psr.append(line.split()[1])
            ra_all_psr.append(float(line.split()[2]))
            dec_all_psr.append(float(line.split()[3]))
            if (line.split()[4] == "-10.0") or (line.split()[4] == "*"):
                p0_all_psr.append(0.0)
            else:
                p0_all_psr.append(float(line.split()[4]))
            # Check for good values of p-dot
            if (line.split()[5] == "-10.0") or (line.split()[5] == "*"): 
                pd_all_psr.append(0.0)
            else:
                pd_all_psr.append(float(line.split()[5]))
            # And good values of DM (also cut out later)
            if (line.split()[6] == "-10.0") or (line.split()[6] == "*"): 
                dm_all_psr.append(0.0)
            else:
                dm_all_psr.append(float(line.split()[6]))
            if len(line.split()) > 7:
                stat_all_psr.append(line.split()[7])
            else:
                stat_all_psr.append('*')
    num_all_psr = np.array(num_all_psr)
    name_all_psr = np.array(name_all_psr)
    ra_all_psr = np.array(ra_all_psr)
    dec_all_psr = np.array(dec_all_psr)
    p0_all_psr = np.array(p0_all_psr)
    pd_all_psr = np.array(pd_all_psr)
    dm_all_psr = np.array(dm_all_psr)
    stat_all_psr = np.array(stat_all_psr)


    north = dec_all_psr > -40 - ang_max
    print "%d pulsars in survey area. Finding closest beams..." %len(num_all_psr[north])
    beam_all_psr = []
    good_beams = np.array([True for n in name_all_psr[north]])
    pulsars=Table.read(file_psr,format='ascii',names=('num','Jname','RAdeg','Decdeg','P0','P1','DM'))
    p_coord=SkyCoord(pulsars['RAdeg'],pulsars['Decdeg'],unit=('deg','deg'))
    for i in range(len(p_coord)):
	beam_all_psr.append(t[p_coord[i].separation(SkyCoord(t['RAdeg'],t['Decdeg']))<0.5*u.deg]['pointing'])



    for i in range(len(beam_all_psr)):
            if len(beam_all_psr[i])==0:
	        good_beams[i] = False
	        continue
	    for j in range(len(beam_all_psr[i])):
	    	beam_all_psr[i][j]=beam_all_psr[i][j].strip('GBNCC')

#Get T_sky information from Joe's file
sky_file = '/users/amcewen/GBNCC_work/skytemp.dat'
read_file2 = open(sky_file,'r')

temp_sky = []
ra_temp = []
dec_temp = []
name_temp = []

for line in read_file2:
    temp_sky.append(float(line.split()[0]))
    ra_temp.append(float(line.split()[1]))
    dec_temp.append(float(line.split()[2]))
    name_temp.append(line.split()[3])

temp_sky = np.array(temp_sky)
ra_temp = np.array(ra_temp)
dec_temp = np.array(dec_temp)
name_temp = np.array(name_temp)
temp_sys = np.array(temp_sky + 70)   # system temperature at 350MHz

read_file2.close()

file_flux = '/users/rspiewak/pulsars/GBNCC_ATNF_fluxes.dat'
psr_flux = np.loadtxt(file_flux,usecols=(1,),dtype=str,unpack=True)
s400_flux,s1400_flux,spindx_flux = np.loadtxt(file_flux,usecols=(10,12,14),unpack=True)
s350_flux = []
for S400,S1400,spindx in zip(s400_flux,s1400_flux,spindx_flux):
    if S400 <= 0 and S1400 <= 0:
        s350_flux.append(0.0)
    elif S400>0 and S1400>0 and spindx<-9:
        s350_flux.append(S400*(0.35/0.4)**(np.log(S400/S1400)/-1.253))
    elif S400>0 and spindx>-9:
        s350_flux.append(S400*(0.35/0.4)**spindx)
    elif S1400>0 and spindx>-9:
        s350_flux.append(S1400*(0.35/1.4)**spindx)
    elif S400>0:
        s350_flux.append(S400*(0.35/0.4)**-1.7)
    elif S1400>0:
        s350_flux.append(S1400*(0.35/1.4)**-1.7)
    else:
        s350_flux.append(0.0)
        
s350_flux = np.array(s350_flux)



name_found = []
name_detect = []
name_unobs = []
ra_found = []
dec_found = []
dist_found = []
snr_e_found = []
snr_found = []
S_found = []
dS_found = []
beam_found = []
MJD_found = []
ra_detect = []
dec_detect = []
beam_detect = []
ra_unobs = []
dec_unobs = []
dist_unobs = []
beam_unobs = []
name_short = []

if proc_psr == True:
    print "Checking S/N for profiles, with threshold = %.1f" %snr_min

beam_rejects = np.array([13691,17237,19064,72172,80142,83626,114632,115242,120582])

# Find observations and, optionally, process with PRESTO
# For single object use:
#for beam in beam_num[lim]:
# For multi-object use:
for beams_psr, name_psr in zip(beam_all_psr,name_all_psr[north][good_beams]):
    print ""
    print ""
    print name_psr, beams_psr
    print ""
    if 'proc' in args:
	proc_psr = True   # Check and correct processing state
    MJD_psr = 10000    # random number
    if 'file' in args:   # Multi-object use
	params = (name_all_psr[north] == name_psr)
	num_psr = num_all_psr[north][params][0]
	ra_psr = ra_all_psr[north][params][0]
	dec_psr = dec_all_psr[north][params][0]
	p0_psr = p0_all_psr[north][params][0]
	p0_short = int(p0_psr*100000 + 0.499)/100.
	pd_psr = pd_all_psr[north][params][0]
	dm_psr = dm_all_psr[north][params][0]
	stat_psr = stat_all_psr[name_all_psr == name_psr][0]
    if name_psr[0] == "C" or stat_psr == 'L' or stat_psr == 'G':
	unpub_psr = True
    else:
	unpub_psr = False
    if f_pub == True and unpub_psr == True:
	continue     # Skip unpublished pulsars when flag used
    if name_psr in name_temp:
	temp_psr = temp_sys[name_temp == name_psr]
    else:
	if ra_psr < 10:
	    lim_t = np.logical_or(ra_temp > 350+ra_psr, ra_temp < ra_psr+10)
	elif ra_psr > 350:
	    lim_t = np.logical_or(ra_temp > ra_psr-10, ra_temp < ra_psr-350)
	else:
	    lim_t = np.logical_and(ra_psr > ra_temp-10,ra_psr < ra_temp+10)
	dist = np.array([ang_offset(ra_psr,dec_psr,RA,DEC) for RA,DEC
			 in zip(ra_temp[lim_t],dec_temp[lim_t])])
	dist,dist_temp = zip(*sorted(zip(dist,temp_sys[lim_t])))
	dist = np.array(dist)
	dist_temp = np.array(dist_temp)
	if len(dist) == 0:
	    temp_psr = 100
	else:
	    for i in range(3):
		if dist[i] == 0:
		    dist[i] = 1e-6
	    temp_psr = np.average(dist_temp[:3],weights=1/dist[:3])
    if dm_psr <= 0 or p0_psr <= 0:   # Check the DM and P for the current pulsar
	proc_psr = False
	print "PSR %s rejected due to invalid parameters" %name_psr
    par_check = False
    par_psr = glob('%sparfiles/%s_short.par' %(home_dir,name_psr))
    if len(par_psr) == 0 and proc_psr == True and len(name_psr)>=10:
	print "Cannot find par file in %sparfiles/" %home_dir
    elif len(par_psr) > 0:
	par_check = True
    if proc_psr == True:
	if len(glob("%s_temp" %name_psr)) == 0 and pull_new == True:
	    os.mkdir("%s_temp" %name_psr)
	elif len(glob("%s_temp" %name_psr)) == 0:
	    print "Skipping new pulsar: "+name_psr
	    continue
	os.chdir("%s_temp/" %name_psr)
	if len(glob('%s_short.par' %name_psr)) == 0 and par_check == True:
	    sproc.call('cp -s %s .' %par_psr[0],shell=True)
    for beam_s in beams_psr:
		beam=int(beam_s)
		if beam in beam_rejects:
                    print "Beam #%d rejected due to known issues with file/obs." %beam
 		    continue
		print beam
		if (dm_psr > 0 and p0_psr > 0) and '-proc' in args:
		    proc_psr = True
		blob = glob('%s*/gu*G????%s_*fits' %(data_dir,beam_s))
		if len(blob) > 0:
#		    rfi_tar = 
		    ra_cand = ra_full[beam]
		    dec_cand = dec_full[beam]
		    MJD_beam = int(blob[0].split('/')[-1].split('_')[1])
		    ra_found.append(ra_cand)
		    dec_found.append(dec_cand)
		    name_found.append(name_psr)
		    dist_found.append(ang_offset(ra_cand,dec_cand,ra_psr,dec_psr))
		    beam_found.append(beam)
		    MJD_found.append(MJD_beam)
		    # Process files to find detections
		    if proc_psr == True:
			if os.getcwd() != work_dir+'%s_temp' %name_psr:
			    if len(glob(work_dir+'%s_temp' %name_psr)) == 0:
				os.chdir(work_dir)
				os.mkdir(name_psr+'_temp')
			    os.chdir(work_dir+'%s_temp/' %name_psr)
			blob = glob('%s2*/gu*G????%s_*fits' %(data_dir,beam_s))
			if len(glob('gu*G????%s_*fits' %beam_s))==0 and len(blob)>0 and pull_new==True:
			    sproc.call("cp -s %s ." %blob[-1], shell=True)
			elif len(glob('gu*G????%s_*fits' %beam_s))==0 and len(blob)>0:
			    print "Skipping new data"
			    continue
			A = sproc.Popen(["ls gu*G????%s_*fits" %beam_s], stdout=sproc.PIPE, shell=True)
			file_beam = A.communicate()[0]
			scan_beam = int(file_beam.split('_')[3])
			rfi_std = "-time 1 -timesig 2 -freqsig 3"
			if "2bit" in file_beam and rfi_fil == True:
			    raw_opt = " -noscales -noweights"
			    rfi_std = "-time 1 -freqsig 3"
			else:
			    raw_opt = " "
			if MJD_beam > 56710 and rfi_fil == True:
			    rfi_opt = " -zapchan 2470:3270"
			else:
			    rfi_opt = " "
			if (len(glob('*_%s_rfifind.mask' %beam_s)) == 0) or \
			   (len(glob('*%s*.bestprof' %beam_s)) == 0 and
			    len(glob('*%s_*.stats' %beam_s)) == 0):
			    rfi_other = glob(work_dir+'*/*_%s_rfifind.stats' %beam_s)
			    for blah in rfi_other:
				if blah.split('/')[-2] != name_psr+'_temp':
				    if len(glob("%s_%s_rfifind.mask" %(name_psr,beam_s))) > 0:
					# Avoid having multiple mask files 
					sproc.call('rm %s_%s*.mask' %(name_psr, beam_s), shell=True)  
				    # Save time on rfifind 
				    sproc.call('cp -s %s.* .' %blah.split('.')[0], shell=True) 
				    rfi_othertxt = glob(work_dir+'*/rfifi*%s_*txt' %beam_s)
				    if len(rfi_othertxt) > 0:
					sproc.call('cp -s %s .' %rfi_othertxt[0], shell=True)
				    break
			    if len(glob('*_%s_rfifind.stats' %beam_s)) == 0:
				print "Running rfifind for %s..." %name_psr
				rfi_out = open('rfifind_%s_output.txt' %beam_s,'w')
				rfi_out.write("nice rfifind %s%s%s -o %s_%s %s"
					      %(rfi_std,raw_opt,rfi_opt,name_psr,beam_s,file_beam))
				p = sproc.Popen("nice rfifind %s%s%s -o %s_%s %s"
						%(rfi_std,raw_opt,rfi_opt,name_psr,beam_s,file_beam),
						stdout=rfi_out, shell=True)
				p.wait()
				rfi_out.flush()
				rfi_out.close()
				del p
			if len(glob('rfifind_%s_output.txt' %beam_s)) > 0:
			    B = sproc.Popen(["grep 'good' rfifind_%s_output.txt | awk '{ print $6 }'"
					     %beam_s],stdout=sproc.PIPE,shell=True).communicate()[0]
			    if B == '' and MJD_beam < 56710:
				bandwidth = 90.0
			    elif B == "" and MJD_beam > 56710:
				bandwidth = 75.0
			    elif len(B) < 4:
				proc_psr = False
				bandwidth = 0.0
			    else:
				bandwidth = float(B[1:-3])
			else:
			    bandwidth = 90.0
			flag_search = '-nosearch'
			if len(name_psr) < 10 and p0_psr < 0.1:
			    flag_search = '-nodmsearch'
			nintx = 64
			if p0_psr > 0.5:
			    nintx = 16
			elif (p0_psr > 0.1) or (p0_psr > 0.5):
			    nintx = 32
			nbinx = nbin
			if p0_psr < 1.7e-3 and nbin > 18:
			    nbinx = 16
			elif p0_psr < 2.5e-3 and nbin > 20:
			    nbinx = 20
			elif p0_psr < 3.5e-3 and nbin > 30:
			    nbinx = 30
			prep_str = "-n %d -nsub %d -npart %d -fine%s -mask *%s*.mask -noxwin %s" \
				   %(nbinx, nsub, nintx, raw_opt, beam_s, file_beam)
			if len(glob('guppi*%s*.pfd.bestprof' %beam_s)) == 0 and proc_psr == True and\
			   unpub_psr == False:
			    fold_out = open('prepfold_%s_output.txt' %beam_s,'w')
			    fold_out.write("prepfold -psr %s %s %s"%(name_psr[1:],flag_search,prep_str))
			    p = sproc.Popen("prepfold -psr %s %s %s" %(name_psr[1:],flag_search,
							    prep_str),shell=True,stdout=fold_out)
			    p.wait()
			    fold_out.flush()
			    fold_out.close()
			    del p
			prof_cand = glob('guppi*%s*.pfd.bestprof' %beam_s)
			if len(prof_cand) == 0 and proc_psr == True and par_check == True:
			    print "Attempting to fold %s with par file..." %name_psr
			    if len(glob("prepfold_%s_output.txt"%beam_s)) > 0:
				sproc.call("rm prepfold_%s_output.txt"%beam_s, shell=True)
			    fold_out = open('prepfold_%s_output.txt' %beam_s,'w')
			    fold_out.write("prepfold -par *t.par %s %s"%(flag_search,prep_str))
			    p = sproc.Popen("prepfold -par *t.par %s %s"%(flag_search,prep_str), 
					    shell=True, stdout=fold_out)
			    p.wait()
			    fold_out.flush()
			    fold_out.close()
			    del p
			    prof_cand = glob("guppi*%s*.bestprof" %beam_s)
			if len(prof_cand) == 0 and proc_psr == True:
			    print "Attempting to fold %s with simple parameters..." %name_psr
			    if p0_psr < 0.1 or len(name_psr) < 10:
				flag_search = "-nodmsearch"
			    if len(glob("prepfold_%s_output.txt"%beam_s)) > 0:
				sproc.call("rm prepfold_%s_output.txt"%beam_s, shell=True)
			    fold_out = open('prepfold_%s_output.txt' %beam_s,'w')
			    fold_out.write("prepfold -p %.11f -pd %.2e -dm %.3f %s %s" 
					   %(p0_psr, pd_psr, dm_psr, flag_search, prep_str))
			    p = sproc.Popen("prepfold -p %.11f -pd %.2e -dm %.3f %s %s" 
					    %(p0_psr, pd_psr, dm_psr, flag_search, prep_str), 
					    shell=True,stdout=fold_out)
			    p.wait()
			    fold_out.flush()
			    fold_out.close()
			    del p
			    prof_cand = glob('guppi*%s*.pfd.bestprof' %beam_s)
			if len(prof_cand) == 0 or proc_psr == False:
			    print "Not fully processed, skipping..."
			    snr_e_found.append(-10.0)
			    snr_found.append(-10.0)
			    S_found.append(-10.0)
			    dS_found.append(10.0)
			else:
			    # Cut extra RFI "manually"
			    if int(file_beam.split("_")[1]) > 56710 and rfi_fil == True \
			       and len(glob('prepfold_%s_*txt')) > 0:
				p = sproc.Popen("tail -4 prepfold_%s_*txt | head -1" %beam_s,
						shell=True,stdout=sproc.PIPE).communicate()[0]
				if float(p.split()[-1]) > 1.3:
				    print "Removing extra RFI"
				    pfd_out = open("show_pfd_%s_output.txt" %beam_s,'w')
				    p = sproc.Popen('show_pfd -killsubs 63,77:102 -noxwin %s*pfd'
						%file_beam.split('.')[0],shell=True,stdout=pfd_out)
				    p.wait()
				    pfd_out.flush()
				    pfd_out.close()
				    del p
			    # Deal with multiple bestprof files
			    blob = np.array(glob('gu*%s*prof' %beam_s)) 
			    if len(blob) == 1:
				prof_file = blob[0]
			    elif len(blob) > 1:
				chisq = np.array([float(sproc.Popen(["grep 'Reduced' %s | \
	awk '{print $5}'" %prof_file], stdout=sproc.PIPE, shell=True).communicate()[0])
						  for prof_file in blob])
				prof_file = blob[chisq == chisq.max()][0]
			    elif len(blob) == 0:
				sys.exit('No bestprof file found for %s' %name_psr)
			    num, val = np.loadtxt(prof_file, dtype='float', unpack=True)
			    lim = val < val.mean()
			    i = 0
			    while i < 3:
				lim2 = val < val[lim].mean()+2.5*val[lim].std()
				lim = val < val[lim2].mean()+2.5*val[lim2].std()
				if val[lim].std() == val[lim2].std():
				    break
				i = i+1
			    W = len(val)-len(val[lim])+0.0
			    if W == 0.0:
				W = 1.0
			    elif W > 1:
				high_val = num[val == val.max()][0]
				for i in range(len(lim)):
				    if i != high_val and lim[i] == False and \
				       lim[i-1] == True and lim[(i+1)%len(lim)] == True:
					lim[i] = True
				W = len(val)-len(val[lim])+0.0
				
			    P = len(val)+0.0
			    sig = val[lim].std()
			    ang_sep = ang_offset(ra_cand,dec_cand,ra_psr,dec_psr)
			    lim_flux = psr_flux == name_psr
			    if (not name_psr in psr_flux) or s350_flux[lim_flux][0]<=0:
				snr_exp = -10.0
			    else:
				snr_exp = s350_flux[lim_flux][0]*2*np.sqrt(2*120*bandwidth)*\
					  np.sqrt((P-W)/W)/(temp_psr*np.exp((ang_sep/0.3)**2 / 1.5))
			    snr_e_found.append(snr_exp)
			    snr_beam = np.abs(np.array([(n-val[lim].mean()) for n in
							val]).sum())/(sig*np.sqrt(W))
			    S_beam = snr_beam*temp_psr*np.sqrt(W/(P-W))/(2*np.sqrt(2*120*bandwidth))
			    S_offset = np.exp((ang_sep/0.3)**2 / 1.5)*S_beam
			    # estimate error in flux using dW=1, dSNR=25%, dT_sys=10K,
			    #  dT_obs=5s, dBW=10MHz
			    # includes estimate of error in angular separation (offset
			    #  from center of beam)
			    dS_offset = np.sqrt((S_beam*np.sqrt((1/4.)**2 + (10/temp_psr)**2 +
						    0.5*(1/24.)**2 + 0.5*(10./bandwidth)**2 +
						    0.5*(P/(W*(P-W)))**2))**2 +
					     (0.1*ang_sep*np.exp((ang_sep/0.3)**2/1.5))**2)
			    if snr_beam < 0.001:
				S_offset = 0.0
				dS_offset = 0.0
			    elif snr_beam < snr_min and snr_beam >= 0.001 and \
			       ang_offset(ra_cand,dec_cand,ra_psr,dec_psr) >= 0.79:
				S_offset = 100.0
				dS_offset = 100.0
			    #chi2_beam = np.array([((n-val.mean())/val.std())**2 for n in val]).sum()
			    snr_found.append(snr_beam)
			    S_found.append(S_offset)
			    dS_found.append(dS_offset)
			    #if chi2_beam > chi2_min:
			    if snr_beam > snr_min:
				print "PSR %s detected in beam #%d on MJD %d with S/N of %.3f; expected S/N %.3f" \
				    %(name_psr, beam, MJD_beam, snr_beam, snr_exp)
				ra_detect.append(ra_cand)
				dec_detect.append(dec_cand)
				name_detect.append(name_psr)
				beam_detect.append(beam)
				if len(glob(work_dir+'detection_plots')) == 0:
				    os.chdir(work_dir)
				    sproc.call('mkdir detection_plots',shell=True)
				    os.chdir(name_psr+"_temp/")
				pfd_str = "detection_plots/guppi_%d_GBNCC%s_%s.pfd" \
					  %(MJD_beam,beam_s,name_psr)
				if len(glob('../%s.*' %pfd_str)) == 0:
				    if len(glob('g*%s*2bit*ps' %beam_s)) > 0:
					file_str = "%s*2bit" %beam_s
				    else:
					file_str = beam_s
				    if len(glob('g*%s*ps'%file_str)) > 1:
					blob = glob('g*%s*prof'%file_str)
					check_f = ''
					check_n = 0.0
					for blah in blob:
					    try:
						cs = float(sproc.Popen(["grep 'Reduced' %s | \
	awk '{print $5}'" %prof_file], stdout=sproc.PIPE, shell=True).communicate()[0])
						if cs > check_n:
						    check_f = blah
					    except:
						pass
					if check_f != '':
					    file_str = check_f[3:-10]
				    os.chdir('%sdetection_plots' %work_dir)
				    sproc.call('cp -s %s%s*/g*%s*ps %s.ps'%(work_dir,name_psr,
						file_str,pfd_str.split('/')[1]), shell=True)
				    sproc.call('cp -s %s%s*/g*%s*prof %s.bestprof' %(work_dir,
						name_psr,file_str,pfd_str.split('/')[1]),shell=True)
				    os.chdir("%s%s_temp" %(work_dir,name_psr))
				if use_comp == 'zuul' and \
				       len(glob("%s%s*" %(home_dir,pfd_str))) == 0:
				    os.chdir("%sdetection_plots/" %home_dir)
				    sproc.call("cp -s %s%s* ." %(work_dir,pfd_str), shell=True)
				    os.chdir("%s%s_temp" %(work_dir,name_psr))
				elif use_comp == 'GB' and len(glob("%s%s*" %(home_dir,pfd_str))) == 0:
				    sproc.call("cp gu*%s_%s*ps %sdetection_plots/"
					       %(beam_s,name_psr,home_dir), shell=True)
				    sproc.call("cp gu*%s_%s*prof %sdetection_plots/"
					       %(beam_s,name_psr,home_dir), shell=True)
			    else:
				if snr_exp <= 0:
				    exptd = "Unknown exp. : %d" %snr_exp
				elif snr_exp > snr_min*1.5:
				    exptd = "Unexpected"
				else:
				    exptd = "Expected"
				print "PSR %s not detected in beam #%d; S/N < %.1f: %s" \
				    %(name_psr, beam, snr_min, exptd)
			os.chdir(work_dir)
		    else:
			snr_e_found.append(-1.0)
			snr_found.append(-1.0)
			S_found.append(-1.0)
			dS_found.append(1.0)
		else:
		    print "Cannot find beam #%d on %s for PSR %s" %(beam, use_comp, name_psr)
		    ra_unobs.append(ra_full[beam])
		    dec_unobs.append(dec_full[beam])
		    name_unobs.append(name_psr)
		    dist_unobs.append(ang_offset(ra_full[beam],dec_full[beam],ra_psr,dec_psr))
		    beam_unobs.append(beam)
if len(name_unobs) > 0 and name_unobs[-1] == name_psr and (len(name_found) == 0 or name_found[-1] != name_psr):
	print "Found no beams on %s for PSR %s" %(use_comp, name_psr)
print ""        
print "YYYEEEEEAAAAAAA (code ran successfully)"
print ""
ra_found = np.array(ra_found)
dec_found = np.array(dec_found)
name_found = np.array(name_found)
dist_found = np.array(dist_found)
snr_e_found = np.array(snr_e_found)
snr_found = np.array(snr_found)
S_found = np.array(S_found)
dS_found = np.array(dS_found)
beam_found = np.array(beam_found)
MJD_found = np.array(MJD_found)
ra_detect = np.array(ra_detect)
dec_detect = np.array(dec_detect)
name_detect = np.array(name_detect)
beam_detect = np.array(beam_detect)
ra_unobs = np.array(ra_unobs)
dec_unobs = np.array(dec_unobs)
name_unobs = np.array(name_unobs)
dist_unobs = np.array(dist_unobs)
beam_unobs = np.array(beam_unobs)

name_short = []
for NAME in name_detect:
    if (len(name_short) > 0 and name_short[-1] != NAME) or len(name_short) == 0:
        name_short.append(NAME)

name_short = np.array(name_short)


## Summarize results ##
print "#####  Summary  #####"
print "%d pulsars in survey area\nDetected %d pulsars in %d beams\nUnable to find \
observations on %s for %d beams" \
    %(len(num_all_psr[north]),len(name_short),len(ra_detect),use_comp,len(ra_unobs))

# Write out summary to file in readable form
jj = max([int(blah.split('_')[-1].split('.')[0]) for blah in
          glob('%sGBNCC_ATNF_beamcheck_%s*.txt' %(home_dir,use_comp))])+1
writefile = open("%sGBNCC_ATNF_beamcheck_%s_%d.txt" %(output_loc,use_comp,jj), 'w')
writefile.write("# Using S/N > %s for detection threshold on %s\n" %(snr_min,use_comp))
writefile.write("# Using beams within %.1f deg of pulsars from %s\n" %(ang_max,file_psr))
if '-pub' in args:
    writefile.write('# Only searched for published pulsars\n')
writefile.write("# Flux adjusted for offset from beam center\n")
writefile.write("# Expected S/N estimated from published flux values\n")
writefile.write("# Flux estimated from measured S/N, uncertainty estimated from dependencies\n")
writefile.write("# File generated %s (ET)\n" %strftime("%Y-%m-%d %H:%M:%S"))
writefile.write("#   PSRJ  \t  Distance  \t   MJD   \t    S/N exp. \t   S/N meas. \t Flux \
Estimate\t Flux Error \t  Status \n")
writefile.write("#\t\t    (deg)   \t  (day)\t\t\t\t\t\t     (mJy)  \t   (mJy)  \n")
for NAME,DIST,MJD,SNR_E,SNR,FLUX,dFLUX,BEAM in zip(name_found,dist_found,MJD_found,snr_e_found,
                                                   snr_found,S_found,dS_found,beam_found):
    writefile.write("%10s\t %.7f \t  %5d  \t %7.3f\t %7.3f\t   %7.3f    \t %7.3f    \t"
                    %(NAME,DIST,MJD,SNR_E,SNR,FLUX,dFLUX))
    if len(beam_detect) > 0 and len(beam_detect[np.logical_and(beam_detect == BEAM,
                                                               name_detect == NAME)]) > 0:
        writefile.write("Detected in Beam #%6d\n" %BEAM)
    else:
        writefile.write("Not Detected in #%6d\n" %BEAM)

writefile.write("# Pulsars with beams not found on %s\n#   PSRJ  \tDistance (deg)\n" %use_comp)
for NAME,DIST,BEAM in zip(name_unobs,dist_unobs,beam_unobs):
    writefile.write("%10s\t  %.7f  \tBeam #%6d Not Found\n" %(NAME,DIST,BEAM))

writefile.close()
#sproc.call('cp -s %sGBNCC_ATNF_beamcheck_%s_%d.txt .' %(output_loc,use_comp,jj),shell=True)


    
"""
# Make plot of area and beams
b = np.linspace(0,2*np.pi,400)
x_circle = ra_psr + 0.6*np.cos(b)
y_circle = dec_psr + 0.6*np.sin(b)

plt.figure(num=1)
fig = plt.gcf()
fig.set_size_inches(4,4)

if len(ra_found) > 0:
    plt.scatter(ra_found,dec_found,c='b',marker=u'o')
plt.scatter(ra_psr,dec_psr,c='g',marker=u'o')
if len(ra_detect) > 0:
    plt.scatter(ra_detect,dec_detect,c='g',marker=u'o')
if len(ra_unobs) > 0:
    plt.scatter(ra_unobs,dec_unobs,c='r',marker=u'o')
plt.plot(x_circle,y_circle,'k--')
plt.xlim(ra_psr-0.8,ra_psr+0.8)
plt.ylim(dec_psr-0.8,dec_psr+0.8)

ax = plt.gca()
ax.invert_xaxis()
plt.xlabel('RA (degrees)')
plt.ylabel('Dec (degrees)')


plt.savefig('beamcheck_%s.png' %name_psr, bbox_inches='tight')
plt.savefig('beamcheck_%s.ps' %name_psr)
"""



