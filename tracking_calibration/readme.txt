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