MATLAB code for calculation of PCOV

1.CalculatingPSDwithMaskinparts (main code)
 Input needed: file path where 2048 IRM frames are saved
               No of cells to be analyzed in one frame
               draw cell boudary, background
               Conversion factor
               Frame rate            
2.psdIMcovMaskpixelreadJG1024Long (function file)

______________
MATLAB code for calculation of temporal, spatial fluctuation, tension, active temp, 
cytoplasmic viscosity, Confinement of all FBRs

1.Multifol_findFBRnew (main code)
Input needed: file path where 2048 IRM frames are saved
               No of cells to be analysed
               Draw a rectangular box around a particular cell to be analyzed
               draw cell boudary, nucleus
               Conversion factor
               Imin
               FBR size
               Frame rate
2.fitPSD (function file)
________________
MATLAB code for calculation of probability distribution of temporal fluctuation and tension

1.consolidate (main code)
Input needed: file path where .mat files generated from "Multifol_findFBRnew" MATLAB code are saved
               
2.compareandplotsd (main code)
 Input needed: load the consolidated .mat file generated from "consolidate" code for generation of probability distribution of temporal fluctuation
3.compareandplot (main code)
Input needed: load the consolidated .mat file generated from "consolidate" code for generation of probability distribution of tension.

_________________
MATLAB code for correlation analysis between Piezo Fluorescence and tension value
 
1. Piezo_tensionloop (main code)
Input needed: give the path where the fluorescence images and tension .mat file(generated from tension mapping) are saved.

__________________
MATLAB code for Tension Mapping

1.tensionmappingsinglefile (main code)
 Input needed: load background pcov file "allpsd3_back_1.mat" and cell file "allpsd3_cell_1.mat" in this format. These files are generated from "CalculatingPSDwithMaskinparts" code for PCOV analysis

2.psdsimfunc6 (function file)
3.calculateR2 (function file)



FBR Size :M=4;N=4;
Frame rate:=19.91
Inmin=11005.22;
conv=284.90;

DATA Set is given for testing of the MATLAB Codes. This single dataset does not include all the data  used for all analysis.