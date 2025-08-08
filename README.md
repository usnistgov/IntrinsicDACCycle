# IntrinsicDACCycle

This reopository introduces tools for assessing the performance of sorbent materials for Direct Air Capture (DAC) from a thermodynamic treatment of the refresh cycle and the intrinsic material properties. 

These tools provide the theoretical upper limit of the CO2 captured per energy, and theoretical upper limit on the purity of the captured CO2.

There are also tools to find Pareto optimum refresh cycles for the sorbents. 

## Contents
```src``` contains the main Julia scripts for perfoming the Intrinsic DAC analysis.
```ScriptsForRunningICD``` contains Julia scripts (as Pluto notebooks) that were run on the HPC environments that performed the Intrinsic DAC analysis.
```FEASSTscripts``` contains the Python scripts that were used to run the FEASST calculations of adsorption properties in the HPC environments.
```Example_FEASST_Scripts``` contains self-contained example scripts for performing FEASST calculations, assuming compatible installations of FEASST.
    - Scripts in ```CO2_Capacity``` require FEASST version 0.22.0
    - Scripts in ```Kh``` require FEAST version 0.19.0
    - Scripts in ```N2_Atmospheric``` require FEASST 0.19.0


## Development

This code base is under active development.

## Coorespondence

Austin McDannald \
austin.mcdannald@nist.gov \
National Institute of Standards and Technology \
Material Measurement Laboratory \
Materials Measurment Science Division \
Data and AI-Driven Materials Science Group 
