# SingleDoorFreezer


#python3.8 #CoolProp /n
This project is a SingleDoor refrigerator device dynamic simulator. It's main goal is to predict the temperautre varia 
tions during device operation./n It "mimics" dynamic simulation by use of steadystate equations for small timesteps and updating key values./n 
In main file(Freezer.py) all device related parameters should be defined. Postprocessing.py is used to read saved results and plot graphs./n
As a base of thermodynamic properties CoolProp is used./n

Task to do:/n
-Implement model of capillary throttling with partial differential equation/n
-Implement moisture handling and defrost algorithm/n
-Implement fram gasket heater model/n
-Improve heat exchanger calculations by division to small pieces/n 
-Optimize using of CoolProp/n
-Improve system "self-regulation"/n
-clean mess in code/n
