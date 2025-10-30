function [ModelP, GenP, EpiP, inputData, X_0] = loadPreprocessedData(inputModelFile, inputDataFile, simStartDate, verbose)
%% 1. INPUTS
    
    % Load Model Parameters
    [ModelP, EpiP, GenP] = loadModelParam(inputModelFile);
    
    %% 3. TRANSFORM LOADED DATA 

    % Number of total population by age group
    GenP.N_pop_AG = [0.21 0.77 0.02] * GenP.Npop;
    GenP.inputFlag = 1;   % To load or not load data
    assert(sum(GenP.N_pop_AG) == GenP.Npop, "Population split into groups is different than the total population");
        
    % 3.1 Load information for the variants modeled
    GenP.nAgeGroups    = length(GenP.N_pop_AG);
    variantNames       = ModelP.infoVar.Variant;          % Variant names
    EpiP.variantNames  = replace(variantNames, "_", "."); % Use proper variant names for legend
    EpiP.nVariants     = length(variantNames);               % Number of variants
    EpiP.variantFactor = ModelP.infoVar.Infectiousness';  % Infectiousness factor per variant
    EpiP.perc_Rec_Var  = ModelP.infoVar.perc_Rec_Var';      % # of people infected and recovered for each variant
    
    % Number of daily contacts (when needed)
    GenP.nContacts = ModelP.nContacts;
    GenP.nContactsM = table2array(GenP.nContacts(:,2:end));
    
    % Immunity protection per variant
    EpiP.pImm = ModelP.pImm;
    EpiP.Imm_fct = table2array(EpiP.pImm(:,2:end));
    
    % 3.2 Epidemiological parameters
    
    % Create a transition matrix of parameters
    
    % Exposed to Infected
    EpiP.fe_iM = repmat(EpiP.fe_i, EpiP.nVariants, 1); 
    EpiP.fe_rM = 1 - EpiP.fe_iM;
    
    % Infected to Recovered
    EpiP.fi_rM = min(EpiP.fi_r ./ EpiP.variantFactor', 1); 
    EpiP.fi_dM = 1 - EpiP.fi_rM;   
    EpiP.fi_rM
    
    % Recovered to Susceptible
    EpiP.fr_sM = repmat(EpiP.fr_s, EpiP.nVariants, 1); 
    EpiP.fr_iM = 1 - EpiP.fr_sM;
    % 3.3 Initial Data 
    
    % Load number of travelers and cases per day
    inputData = loadData(inputDataFile);
    % Load all cases
    inputData.CasesAdj.Date = days(inputData.CasesAdj.Date- simStartDate);
    rowID = find(inputData.CasesAdj.Date <=0, 1);
    Inf_0 = table2array(inputData.CasesAdj(rowID, 2:end));
    E = GenP.perc_E .* reshape(Inf_0, [GenP.nAgeGroups EpiP.nVariants])';
    I = (1-GenP.perc_E) .* reshape(Inf_0, [GenP.nAgeGroups EpiP.nVariants])';
    
    % Calculate Recovered compartment for each variant per age group
    R = (GenP.N_pop_AG - sum(E+I)) .* EpiP.perc_Rec_Var';
    
    % Calculate initial susceptible people to ALL variants
    S = GenP.N_pop_AG - sum(R) - sum(E+I);
    
    % Define initial deaths
    D = zeros(EpiP.nVariants, length(S));
    
    % Linearize to append vectors
    E_lin = reshape(E', 1, numel(E));
    I_lin = reshape(I', 1, numel(I));
    R_lin = reshape(R', 1, numel(R));
    D_lin = reshape(D', 1, numel(D));
    X_0 = [S E_lin I_lin R_lin D_lin];
    
    % Check that the input sum is equan to the whole population
    assert((sum(X_0) - GenP.Npop) < 1e-3, "Population split into compartments is different than the total population")
    
    % Print message if desired
    if verbose == 1
        fprintf("The model evaluates %i variants and %i age groups.", + ... 
                EpiP.nVariants, GenP.nAgeGroups);
        fprintf("\n-------------------------------------------------\n");
        fprintf("The initial number of susceptible people are:");
        S 
        fprintf("The initial number of infected people per variant \n" + ...
                "and per age group are:");
        E 
        fprintf("The initial number of infected people per variant \n" + ...
                "and per age group are:");
        I 
        fprintf("The initial number of recovered people are:");
        R
    end
    
    % Raw probability of infection by a contact for variant of reference
    EpiP.pInf_V = min(GenP.pInf * EpiP.variantFactor, 1);
    
    % Print message if desired
    if verbose == 1
        fprintf("The infection rate per variant is: ");
        for i = 1:length(EpiP.pInf_V)
            fprintf("Variant#%i infectiousness: %.2f. \n", i, EpiP.pInf_V(i));
        end
    end
    
    % Convert dates into number for better handling during the simulation
    inputData.TravelersDates = days(inputData.Travelers.Date - simStartDate);
    inputData.TravelersM = round(inputData.TravelersM .* GenP.pVis_inf, 0);
    