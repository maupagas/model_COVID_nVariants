% calcRates calculates the rates of various transitions between states (e.g., infection, recovery, death)
%
% Parameters:
%   S: number of susceptible individuals
%   I: number of infected individuals (matrix with age groups in columns)
%   R: number of recovered individuals
%   pinf_V: infection probability for each age group
%   nContactsM: number of contacts per age group
%   EpiP: structure containing parameters of the epidemiological model
%   pImm: proportion of recovered individuals that lose immunity
%   reinf_Flag: flag indicating whether reinfection is allowed
%
% Returns:
%   r: a structure containing the rates of infection, reinfection, recovery, and death
%
% Example:
%   S = 1000; % number of susceptible individuals
%   I = [100 200 300]; % number of infected individuals in each age group
%   R = 100; % number of recovered individuals
%   pinf_V = [0.1 0.2 0.3]; % infection probability for each age group
%   nContactsM = [10 15 20; 5 10 15; 2 5 10]; % number of contacts per age group
%   EpiP.fi_rM = [0.8 0.9 0.95]; % fraction of infected individuals that recover in each age group
%   EpiP.ti_

% function [r] = calcRates(t, S, E, I, R, pinf_V, nContactsM, EpiP, nInfMax, pImm, ...
%     reinf_Flag, InputData)

function [r] = calcRates(t, S, E, I, R, GenP, EpiP, SimP, InputData)
% Fraction of infectious individuals interacting
% per age group (assuming nobody tested and there
% is no reduction in testing)
finf_M = I ./ (S + sum(E + I + R));

% Calculate infection rate
% Infectiousness of infectious individuals interacting with susceptible
% individuals to infection
%     inf_I = bsxfun(@times, sum(nContactsM, 1) .* finf_M, pinf_V');
inf_I = sum(GenP.nContactsM, 1) .* finf_M .* EpiP.pInf_V';

r.r_inf   = inf_I .* S;

% Allow for reinfection or not
if SimP.reinf_Flag

    % Link the immunity to the maximum expected peak of cases
    if GenP.nInfMax == 0
        fImm = 1;
    else
        nInf = sum(sum(E + I));
        if nInf > GenP.nInfMax
            fImm = 0;
        else
            fImm = 1 - nInf/GenP.nInfMax;
%             fImm = fImm_lin/(1+exp(-fImm_lin));
        end
        %             fImm = exp(- nInf/nInfMax);
        %             fImm = (1-exp(-t/365));
        %             nInfV = sum(E+I, 2);
        %             fImm = 0.9995.^nInfV;
    end

    % Calculate the reinfection rates to compartment j (r.r_reinf) and
    % the reinfection rates coming from compartment i (r.r_reinf_0)
    r.r_reinf = EpiP.fr_iM .* inf_I .* (EpiP.Imm_fct .* fImm * R);
    r.r_reinf_src = EpiP.fr_iM .* (fImm' .* EpiP.Imm_fct' * inf_I) .* R;

else
    r.r_reinf = 0;
end

% Transition rate from Exposed to Infected
r.re_i = E .* EpiP.fe_iM ./ EpiP.te_i;
r.re_r = E .* EpiP.fe_rM ./ EpiP.te_i;

% Calculate transition rate from Infected to recovered
r.ri_r = I .* EpiP.fi_rM ./ EpiP.ti_r;
r.ri_d = I .* EpiP.fi_dM ./ EpiP.ti_r;

% Calculate transition rate from Recovered to Susceptible
r.rr_s = R .* EpiP.fr_sM ./ EpiP.tr_s;

% Calculate the number of cases imported
% if isvarname("InputData")
    % Find the rate of arrivals
    rowID = find(InputData.TravelersDates >=floor(t), 1);
    rArr = InputData.TravelersM(rowID, :);
    % Find the rate of departures (equal to rate of arrivals 10 (or tVisit) days before)
    rowID = find(InputData.TravelersDates >=floor(t) - 10, 1);
    rDep = InputData.TravelersM(rowID, :);

    % Reshape to appropriate dimensions
    if isempty(rArr)
        rArr = zeros(size(I));
        rDep = zeros(size(I));
    else
        rArr = reshape(rArr, size(I'))';
        if isempty(rDep)
            rDep = zeros(size(I));
        else
            rDep = reshape(rDep, size(I'))';
        end
    end

    % Store in the r (rates) structure
    r.rArr = rArr;
    r.rDep = rDep;
% end

end