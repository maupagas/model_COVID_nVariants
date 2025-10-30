% Define files to be loaded
inputDataFile = "../input-files/data_AD_model.xlsx";
inputModelFile = "../input-files/modelParameters.xlsx";
SimP.verbose = 1;
SimP.reinf_Flag = 1;
% Define Simulation Parameters
simStartDate = datetime("01/08/2022", "format", "dd/MM/uuuu");
simLastDate  = datetime("01/03/2023", "format", "dd/MM/uuuu");

% Load preprocessed data
[ModelP, GenP, EpiP, inputData, X_0] = ...
    loadPreprocessedData(inputModelFile, inputDataFile, simStartDate, SimP.verbose);

% Prepare and Run Simulation

% Options for the ode solver to limit the maximum step
ode_opts = odeset('MaxStep', 1);
SimP.tFinal = days(simLastDate - simStartDate);
simTime = 0:1:SimP.tFinal;
SimP.simDates = simStartDate + simTime;

 
% Run simulation and time it
tic;
[tsim, X_t] = ode45('modelSimulator', simTime, X_0, ode_opts,...
    GenP, EpiP, SimP, inputData);
toc;

% Check conservation of population over whole simulation
assert(all(round(sum(X_t,2)/1000, 0) == GenP.Npop/1000), "Simulation results of the total populations do not remain constant.")

% Plot results 
plotOutput(X_t, SimP.simDates, inputData, GenP, EpiP)


%% Calculate rates of infection
[~, dailyCases] = cellfun(@(t,x)  modelSimulator(t, x.', ode_opts, GenP, EpiP, SimP, inputData),...
        num2cell(tsim), num2cell(X_t, 2),'uni',0);
% Sum all new infectiosn per at each cell vector 
totalDailyCases = cellfun(@(x) sum(x, 'all'), dailyCases);
