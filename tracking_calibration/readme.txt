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


TRACKING:
For all tracking files and directions, please refer to:

https://github.com/treelinemike/tracker-serial-interface/tree/main

Mike Kokko has created a server and written documentation for collecting Aurora EM tracking data into MATLAB. 

The group has successfully tested the server and collected data with the coil. No data was saved during this test run. 
Please refer to auroura_example.m for the code used during this test run. 
