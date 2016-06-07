(1) Overview:
    Included are utilities for quick meshing of the Cascadia Community Velocity
    Model (CCVM).  The goal of these utilities is to make an okay mesh in a
    short relatively short amount of time.
    
(2) Workflow:

  (a) Unfortunaely I'm pressed for time.  Basically, you have to:
       (i)  Get the latest and greatest CVMs from Art Frankel and name
            them:
              Vpl1c.bin
              Vpl2c.bin
              Vpl3c.bin
              Vsl1c.bin
              Vsl2c.bin
              Vsl3c.bin 
       (ii) Modify the cvm.ini file to encapsulate your region of interest
       (iii) Cross your fingers and run cvm_mesher.  

  (b) Eventually there will be options for converting those wonderful 
      CVM binary files hdf5 so you don't need to know all the hardwired
      constants in include/cvm_constants.h.  Also you could stick an 
      xdmf file next to it and actually look at the model with Paraview.
      cvm_repack would be a good name for that.  cvm_repack'd models could
      then be extracted to produce models for NonLinLoc.  This is because
      if you want to work with real data (why?) you'll probably discover
      the 1D catalog earthquake locations and the 3D earthquake locations
      are different - i.e. you'll have to relocate with a 3D locator like
      NonLinLoc.  cvm2nll would be a good name for that.   And cvm_mesher
      shouldn't do everything I just described because as it stands now
      you have to go in an uncomment things in random spots.

(3) Installation:
  (a) Dependencies
      -  Download and install GSL >= v2
         http://www.gnu.org/software/gsl/
      -  Download and install libiscl:
         https://github.com/bakerb845/libiscl
      -  Download and install zlib
      -  Download and install hdf5
      -  If you want to add topography to your mesh download the following
         netCDF4 topography file:
         https://drive.google.com/folderview?id=0B_12SLGjHpOzSzl6OWh0TUJ0bGc&usp=sharing 
      -  A quick + easy utility for writing vtk files that I used while
         debugging the tripling layer utility is:
         https://github.com/certik/visit_writer
         You additionally will have to enable by setting the 
         -DVISIT_WRITER pre-processing flag. 

