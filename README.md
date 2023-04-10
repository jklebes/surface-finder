# transform_extractsurfaces

The 3-in-1 program transform_extractsurfaces.m handles everything from 3D raw data to selected 2D slices.  The transform, surface finding, and cross section steps have been merged so that 3D data only has to be read into memory once during the analysis.

## transformation

Images from the microscope, taken at a 45deg angle and with a spacing or 4 pixels between planes, need to be transformed via Java/imageJ transformJ affine transform.  Java/imageJ was found to do the affine transform much faster than matlab.

This step can either be done seperately, in which case transformed 3D data and flag ```transfrom=false``` is given.  Or, with ```transfrom=true```,  java commands equivelant to the imageJ script are used to read in, transform, and return the 3D data for further processing.

## surface finding

Surface finding algorithm by Guillermo.

## cross sections

zy slices of the transfromed 3D data are extracted and saved as images for further inspection.  This was combined into the same program so that 3D data doesn't need to get read into emmory twice.
