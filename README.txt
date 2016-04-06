(1) Overview:
    Included are utilities for quick meshing of the Cascadia Community Velocity
    Model (CCVM).  The goal of these utilities is to make an okay mesh in a
    short amount of time.
    
(2) Workflow:

  (a) cvm_resample:

        The workflow usually begins by downsampling and interpolating the CCVM.
        To do this in the ini file set the dx_resample_cvm, dy_resample_cvm, and
        dz_resample_cvm to the desired values (preferably all equal).  Then,
        set the region of interest if modeling a subset of the CVM by specifying
        the lower left corner bound (lat0,lon0) and upper right corner bound 
        (lat1,lon1).  The mesh will be resampled and saved to and hdf5 file.
        An accompanying xdmf file will be produced so you can view the model
        with Paraview.

  (b) cvm_mesher:

        After

  (c) cvm_addtopo:
        With the simple mesh created it can now be deformed to topography.  For
        a topographic model  
         
         
(3) Installation:
  (a) Dependencies
      -  Download and install hdf5.  Note this may require an additional
         library, zlib

