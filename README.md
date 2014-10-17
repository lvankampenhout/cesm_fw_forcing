cesm_fw_forcing
===============

Custom code that implements freshwater (FW) forcing in CESM. Questions to l.vankampenhout[at]uu.nl

User provides freshwater fluxes (liquid, solid) in 12-monthly files. The fluxes are read by POP and inserted (added) at every timestep. Interpolation is applied when appropriate. As the forcing files are given for each calendar year, 3 files are needed in total: the fluxes for the current model year, the fluxes of previous year (needed for interpolation in early January) and the fluxes of next year (needed for interpolation in late December). The simulation needs to be restarted at every new calendar year to update these file pointers.  
Variable names for the liquid and solid fluxes are FW_liquid and FW_solid, respectively. These are hardcoded. The liquid flux is added to POP field ROFF_F and the solid flux to IOFF_F. 

Optionally a blanking can be applied before the freshwater fluxes are added. In the blanking step the fields ROFF_F and IOFF_F are set to zero on a user-defined mask. In this way a flux-replacement design can be achieved. Cells are blanked when mask > 1.0. The variable name for the input file is hardcoded to BLK_MASK.

Limitations
------------

* Only monthly data is fully implemented and tested
* Code must be stopped at every model year (01-Jan) to accomodate file pointer update
* Fixed variable names in forcing and blanking files.

Warnings
--------
* Annual forcing (data_type='annual') has not been fully implemented yet. Amongst others, this only uses the liquid fluxes (variable FW_liquid) and does not apply interpolation. Use it at your own risk!

Compatibility
-------------
These files are based on CESM 1.1.2 but may work with other releases as well.

Files
-----
The following files are present:

* build-namelist
  includes freshwater namelist in the build scheme
  
* namelist_definition_pop2.xml
  defines freshwater namelist (forcing_imau_nml)
  
* forcing_coupled.F90
  implements the file reading and adding of the freshwater fluxes. Optionally the 
  original fluxes are blanked. 

* `forcing_tools.F90`
  adapted interpolation scheme, needed for monthly data. The interpolation for beginning of January
  needs the values from the previous year. Likewise the interpolation of end-December needs the fields
  for the next year. 

How to use
----------
Copy the above files to the SourceMods directory $CASE/SourceMods/src.pop2/

The freshwater forcing must be enabled through the namelist. CESM version 1.1 or higher: 

    cat >> user_nl_pop2 << _EOF
    imau_filename&forcing_imau_nml        = '/home/kampe004/data/fw/FW_GrIS_AIS_1850.nc'
    imau_filename_prev&forcing_imau_nml   = '/home/kampe004/data/fw/FW_GrIS_AIS_1849.nc'
    imau_filename_next&forcing_imau_nml   = '/home/kampe004/data/fw/FW_GrIS_AIS_1851.nc'
    imau_data_type&forcing_imau_nml       = 'monthly'
    imau_blanking&forcing_imau_nml        = .true.
    imau_filename_blanking&forcing_imau_nml   = '/home/kampe004/data/blanking_mask/blk_msk_gx1v6.nc'
    _EOF


