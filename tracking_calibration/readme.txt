CAD FILES:
All necessary CAD files can be found at:
https://cad.onshape.com/documents/2c4d481847ee4e77d08d8edb/w/dc69158cdf2199aaff2963a7/e/71a8c64b8cfca10df6c83b66?renderMode=0&uiState=647cb472180bea11e237f8ad

This linked document includes CAD files for the original probe design and the new probe designed in ENGS 169. The necessary STL files can be found within each probe's respective folder/test-date folder. 

CALIBRATION:
The following files are from calibration protocol for both the new flexable probe and old probe. 
Both calibration protocols are run with the coil being taped to the outside of the probe.

In order to use tracking within either probes, the PCBtoSecondCoil transformation matrix must be updated.
To update the matrix, 
1) Place probes in the center of the aurora target plate. 
2) Collect the SecondCoil matrix and copy into MATLAB.
3) Copy the CalPlatetoEM into MATLAB
4) Calculate the PCBtoSecondCoil matrix by the following relationship: T_plate_to_coil = inv(T_coil_to_em)*T_plate_to_em

For each probe, all the files necessary are included in the respective folder. 

Please refer to Calibration_Protocol.pdf for detailed explanation of how to calibrate a new probe.

Validation:
Once calibrated, provided is a grid with known spacing to determine your registration error.

Print out: calib.io_circles_200x150_8x10_18_8.33.pdf 
Place probe different points.
Save the tracking data for these points and compare with known distances.

Due to printing, the scale might be slightly off. ALWAYS confirm circle diameter and spacing with a ruler. 
It is recommended to ALWAYS use the grid to validate tracking accuracy.


TRACKING:
For all tracking files and directions, please refer to:

https://github.com/treelinemike/tracker-serial-interface/tree/main

Mike Kokko has created a server and written documentation for collecting Aurora EM tracking data into MATLAB. 

The group has successfully tested the server and collected data with the coil. No data was saved during this test run. 
Please refer to auroura_example.m for the code used during this test run. 
