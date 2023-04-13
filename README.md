! This git repo was merged into https://github.com/cjwlab/DSLMDataAnalysis/tree/jason as subfolder, development here is done

# transform_extractsurfaces

The main 3-in-1 program transform_extractsurfaces.m handles everything from 3D raw data to selected 2D slices.  The transform, surface finding, and cross section steps have been merged so that 3D data only has to be read into memory once during the analysis.

## setup

For the ```transform=true``` java route relevant java jars, 

    ij-??.jar
    
    ImageScience.jar
    
must be added to matlab's ```javaclasspath``` (see https://uk.mathworks.com/help/matlab/matlab_external/static-path-of-java-class-path.html) .
You can probably find the jars within the preexisting imageJ/Fiji folder.  

## transformation

Images from the microscope, taken at a 45deg angle and with a spacing or 4 pixels between planes, need to be transformed via Java/imageJ transformJ affine transform.  Java/imageJ was found to do the affine transform much faster than matlab.

This step can either be done seperately, in which case transformed 3D data and flag ```transform=false``` is given.  Or, with ```transform=true```,  java commands equivelant to the imageJ script are used to read in, transform, and return the 3D data for further processing.

The ```transform=true``` Java route suffers from a memory problem.  I deallocate all arrays and Java objects, so that they are in principle flagged for garbage colelction.  Garbage collection does seem to happen some of the time, with amount of memory used on each worker sometimes returning to the original and total memory use approximately holding steady once memory pressure is high.  However, eventually the limit 100% memory usage is hit.  The java garbage collector does not reliably clear everything in time.  The workaround is to restart all workers after a few parfor iterations.

Using the ImageJ script in a seperate preprocessing step has the disadvange that transformed 3D data is written to, read from tif files in the ```_transform_tmp``` directories, using more time and space.
  
## surface finding

Surface finding algorithm by Guillermo used in ```calculate_heightmap``` function.
The image is divided into ```sq_side x sq_side``` (here 25 x 25) square grids ```overlaps``` (here 3) different ways, offset by ```sq_size/overlaps``` (here 8 pixels) from each other.  At each square, the height at which square z-gradient is greatest is chosen.  From among the results from the three different griddings, the image pixels with greatest intensity are chosen.  The final heightmap is smoothed with a box blur.

## cross sections

zy slices of the transfromed 3D data are extracted and saved as images for further inspection.  This was combined into the same program so that 3D data doesn't need to get read into emmory twice.
