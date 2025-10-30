clear all

% Define name for the sensitivity analysis variable
sensitivityInputFileName = "../input-files/modelParameters_sensitivity.xlsx";
inputDataFile = "../input-files/data_AD_model.xlsx";
simStartDate = datetime("01/08/2022", "format", "dd/MM/uuuu");
simLastDate  = datetime("01/03/2023", "format", "dd/MM/uuuu");
SimP.verbose = 0;
SimP.reinf_Flag = 1;

% Save Sensitivity Analysis File Flag
saveSensitivityResults = 1;

% Load preprocessed data
[ModelP, GenP, EpiP, inputData, X_0] = ...
    loadPreprocessedData(sensitivityInputFileName, inputDataFile, simStartDate, SimP.verbose);

% Options for the ode solver to limit the maximum step
ode_opts = odeset('MaxStep', 1);
SimP.tFinal = days(simLastDate - simStartDate);
simTime = 0:1:SimP.tFinal;
SimP.simDates = simStartDate + simTime;


%% Start Montecarlo simulations

% Definite input values to the model to conduct sensitivity analysis on
% infectiousness values

% Load sensitivity analysis file
sensitivity_input_files = strcat('lhs_2/',cellstr(ls("lhs_2/lhs*")));

for id_file=1:length(sensitivity_input_files)

    lhs_file = sensitivity_input_files{id_file};

    output_folder=strcat('output-lhs/', replace(lhs_file, '.mat', ''));
    output_folder=replace(output_folder, 'lhs_2/', '');

    mkdir(output_folder)

    sprintf('Loading %s... \n', lhs_file)

    load(lhs_file)

    infectiousness = X(:,1:5);
    immunity_factors = X(:, 6:end);
    num_immunity_factors = width(immunity_factors);

    numSimulations = length(X);
    variantNames = ModelP.infoVar.Variant;
    numVariants = length(variantNames);

    % Store original values for updating later
    orig_pInf_V = min(GenP.pInf * ModelP.infoVar.Infectiousness', 1);
    orig_Imm_fct = EpiP.Imm_fct;
    % Preallocate output to make interval confidence of results
    simResults = zeros(length(simTime), numVariants, numSimulations);
    totalCases = zeros(numSimulations, 1);

    % Run Montecarlo Simulations for each infectiousness parameter sampled
    for i = 1:numSimulations

        % Update Epidemiological parameter for variant
        EpiP.pInf_V = min(GenP.pInf * infectiousness(i,:), 1);

        % Update immunity factor matrix
        Immunity_matrix = zeros(size(EpiP.Imm_fct));
        immunity_factor_vector = immunity_factors(i,:);
        % Fill the diagonal with the immunity values
        index = 1;
        for k = 2:5
            for j = 1:(k-1)
                if index <= num_immunity_factors
                    Immunity_matrix(k, j) = immunity_factor_vector(index);
                    index = index + 1;
                end
            end
        end

        EpiP.Imm_fct = Immunity_matrix;
        



        % Run Simulation with increased value
        [tsim, X_t] = ode45('modelSimulator', simTime, X_0, ode_opts,...
            GenP, EpiP, SimP, inputData);

        %% Calculate rates of infection
        [~, dailyCases] = cellfun(@(t,x)  modelSimulator(t, x.', ode_opts, GenP, EpiP, SimP, inputData),...
            num2cell(tsim), num2cell(X_t, 2),'uni',0);


        % Sum all new infectiosn per at each cell vector
        totalDailyCases = cellfun(@(x) sum(x, 'all'), dailyCases);
        totalCases(i) = sum(totalDailyCases);

        % Get Total Cases Per Variant
        [totalSimCasesPerVariant] = getTotalCasesPerVariant(X_t, EpiP, GenP);

        % Store Results
        simResults(:,:, i) = totalSimCasesPerVariant;
        % Reset EpiP.pInf_V to their original values
        EpiP.pInf_V = orig_pInf_V;

        % Print counter every 10 simulations

        if mod(i, 250) == 0
            fprintf("Simulation #%i out of %i finished.\n", i, numSimulations);

        end

    end

    % Save results after simulation
    if saveSensitivityResults
        currentDate = yyyymmdd(datetime('today'));
        saveFilename = sprintf("%s/%i_simResults_%s", output_folder, currentDate, replace(lhs_file, 'lhs_2/', ''));
        save(saveFilename)
    end

end