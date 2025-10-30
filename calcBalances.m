% calcBalances calculates the changes in the numbers of susceptible, infected,
% recovered, and deceased individuals over time.
%
% Parameters:
%   r: a structure containing the rates of various transitions between states 
% (e.g., infection, recovery, death)
%
% Returns:
%   dX_dt: a 1x4 vector containing the changes in the numbers of susceptible,
% infected, recovered, and deceased individuals, respectively
%
% Example:
%   r.r_inf = 0.01; % infection rate
%   r.r_reinf = 0.005; % reinfection rate
%   r.ri_r = 0.1; % recovery rate for infected individuals
%   r.ri_d = 0.01; % death rate for infected individuals
%   dX_dt = calcBalances(r); % [dS_dt dI_dt dR_dt dD_dt]

function [dX_dt] = calcBalances(r)
    
    % Balance of all compartments
    dS_dt = -sum(r.r_inf) + sum(r.rr_s, 1);
    dE_dt = +r.r_inf + r.r_reinf - r.re_i - r.re_r;
    dI_dt = +r.re_i  - r.ri_r - r.ri_d + r.rArr;
    dR_dt = +r.ri_r  - r.r_reinf_src - r.rr_s - r.rDep;                 
    dD_dt = +r.ri_d;

    % Linearize to append vectors
    dE_dt_lin = reshape(dE_dt', 1, numel(dE_dt));
    dI_dt_lin = reshape(dI_dt', 1, numel(dI_dt));
    dR_dt_lin = reshape(dR_dt', 1, numel(dR_dt));
    dD_dt_lin = reshape(dD_dt', 1, numel(dD_dt));

%     dXvar_dt = [dI_dt dR_dt dD_dt];
%     dXvar_dt_lin = reshape(dXvar_dt, 1, numel(dXvar_dt));
%     dX_dt = [dS_dt dXvar_dt_lin];
  
    % Merge all vectors together
    dX_dt = [dS_dt dE_dt_lin dI_dt_lin dR_dt_lin dD_dt_lin];

end