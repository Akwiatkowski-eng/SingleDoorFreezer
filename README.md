# SingleDoorFreezer

This project is a SingleDoor refrigerator device dynamic simulator. It's main goal is to predict the temperautre varia 
tions during device operation. It "mimics" dynamic simulation by use of steadystate equations for small timesteps and updating key values. 
In main file(Freezer.py) all device related parameters should be defined. Postprocessing.py is used to read saved results and plot graphs.
As a base of thermodynamic properties CoolProp is used.

Task to do:
-Implement model of capillary throttling with partial differential equation
-Implement moisture handling and defrost algorithm
-Implement fram gasket heater model
-Improve heat exchanger calculations by division to small pieces 
-Optimize using of CoolProp
-Improve system "self-regulation"
-clean mess in code
