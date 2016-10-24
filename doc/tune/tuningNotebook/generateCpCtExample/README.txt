Example of generating a Cp/Ct table for use with FLORIS via a FAST model
Paul Fleming
Aug 11, 2016

A key turbine input required by FLORIS is a Cp/Ct table against wind describing the the turbines.  

This can be derived given a FAST model of a particular turbine, which is demonstrated here.

In this example, a Cp/Ct vs wind are found for the NREL 5-MW reference turbine.

In the subfolder stepSim, the NREL 5MW turbine is set to run in a HH wind file which steps the wind across the range of operation.

Note to re-run this simulation, you need to download FAST (and specifically a compiled version, or compile yourself, where FAST is set to use an external controller)
Used in this work is: FAST_v7.01.00a-bjj_AeroDyn_v13.00.01a-bjj_BladedDLLInterface.exe

The output file needs to include these variables from FAST:
 HorWindV, RotSpeed, BlPitch1, GenTq, RotCp, GenCp, RotCt

Next the ipython notebook calculateCpCtFromFAST.ipynb can be used to produce the tables
The notebook is designed to step through the process such that it could be repeated in other codes if desired

NOTE1: RotCt in FAST7 shouldn't be used directly since it includes gravity loads, the python script accounts for this
       In FAST8 I believe it is possible to directly produce aero Ct

NOTE2: Cp and Ct ultimately produced are "rotor" meaning power on low-speed shaft, not including efficiency losses