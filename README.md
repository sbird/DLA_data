DLA_data
========

DLA data - A collection of datasets for DLAs

Contains:

dla_data.py: Routines for plotting the column density distribution function data

Various CDDF datafiles (citations are inside dla_data.py):
005CDF_z3vwotuv2all.dat
2fn_sdss_dr5.dat       
dndx.txt               
fhix.dat               
fn_celine_z3.dat       
fn_sdss_dr5.dat        
not_2012.dat           
peroux05_z3.dat        
prochaska_05.dat       
prochaska_lls.dat      
summary.dat            


vel_data.py: Routines  for plotting the velocity width and metallicity data, as well as correlations between them

apj469315t2_mrt_mod.txt
apj469315t2_mrt.txt:
Data from Neeleman, M et al 2013: 1303.7239 
velocity widths and metallicities for HIRES spectra. _mod.txt is altered to
be easier for numpy to read.

apj433829t2_mrt.txt    
apj433829t3_mrt.txt    
DLA metallicities from Neeleman, M et al 2012: 1205.5047
Compendium of a lot more metallicities, without reliable velocity widths,
over a wide redshift range. t2 is new to that paper, t3 is literature.
