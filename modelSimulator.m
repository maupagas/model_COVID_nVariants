function [dX_dt, rinf_t] = modelSimulator(t, x, u, GenP, EpiP, SimP, InputData)
  
    % Reshape the input as per the format required
    [S, E, I, R, D] = reshapeInput(x, EpiP.nVariants, GenP.nAgeGroups, 1);

    if SimP.verbose && t == 0
        fprintf("The initial number of susceptible people are:");
        S 
        fprintf("The initial number of exposed people per variant \n" + ...
                "and per age group are:");
        E 
        fprintf("The initial number of infected people per variant \n" + ...
                "and per age group are:");
        I 
        fprintf("The initial number of recovered people are:");
        R
    end

    % Correct any zero-crossing
    if any(S), S(S<0) = 0; end
    if any(E), E(E<0) = 0; end
    if any(I), I(I<0) = 0; end
    if any(R), R(R<0) = 0; end
    
    % Calculate the transition rates for the odes
%     r = calcRates(t, S, E, I, R, pinf_V, nContactsM, EpiP, GenP.nInfMax, Imm_fct, ...
%                 reinf_Flag, data_AD);
    r = calcRates(t, S, E, I, R, GenP, EpiP, SimP, InputData);
    rinf_t = r.r_inf + r.r_reinf;

    if SimP.verbose && t == 0
        for i = 1:EpiP.nVariants
            fprintf("The rate of infection for variant %i per age group are:", i);
            r.r_inf(i,:)
        end
    end
    % Calculate the derivatives
    [dX_dt] = calcBalances(r)';
  
end