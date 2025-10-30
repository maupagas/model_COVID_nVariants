# Model Framework to Simulate _n_ Variants

In this repository, a framework to simulate the spread of an infectious disease that can have multiple variants is presented. 

The model is generalised in such a way that can be segregated into social groups (e.g. age groups, working/non-working class, etc.).

## How to run this model

This framework is distributed in the following manner:

- `input-files`: Folder with all the excels used to run the simulation.
    - `data_AD_model.xlsx`: Contains all the raw input data.
    - `modelParameters.xlsx`: Parameters used to run a single simulation batch.

- `run-files`: Folder with all the MATLAB scripts to run the simulations: 
    - `runSimulationBatch.m`: Runs a simulation with all the parameters defined in `modelParameters.xlsx` from the input-files folder.
    - `runSensitivity_LHS_all.m`: Runs a sensitivity analysis based on LHS.
    

- `output-files`: Stores the results of the sensitivity analysis, the calibration procedure and the figures used in the manuscript in a .fig format for further postprocessing (if required).

## How to run this model

For all runs, the user must add to the MATLAB path the folder `model-COVID-nVariants` and navigate to the `run-files` folder.

### Run a single simulation
For a single run of the model,execute the file `runSimVariants.m`. The parameters of the model can be changed at `./input-files/modelParameters.xlsx` and `data_AD_model.xlsx`.

### Run Sensitivity Analysis
To conduct a sensitivity analysis, execute the file `runSensitivity_LHS_all`. Initial values for sensitivity analysis are stored in the file `./input-files/modelParameters_sensitivity.xlsx`. 

