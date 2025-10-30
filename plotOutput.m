function []=plotOutput(X_t, simDates, InputData, GenP, EpiP)


% Column names to compare
varCols = strcat(strrep(EpiP.variantNames, ".", "_"), "__");
sumVar = zeros(length(InputData.CasesAdj.Date), length(varCols));
for k = 1:length(varCols)
    mask = startsWith(InputData.CasesAdj.Properties.VariableNames, varCols(k));
    Fcols = InputData.CasesAdj(:,mask);
    sumVar(:,k) = sum(table2array(Fcols), 2);
end

colPoints = lines;

% Reshape output to proper format
[S_t, E_t, I_t, R_t, D_t] = reshapeInput(X_t', EpiP.nVariants, GenP.nAgeGroups, 0);

%% Plot totals for all variants
figure(1)
plot(simDates, [sum(S_t, 2), sum(E_t, 2), sum(I_t, 2), sum(R_t, 2), sum(D_t, 2)], 'LineWidth', 1.5)
legendNames = {'Susceptible', 'Total Exposed', 'Total Infected', 'Total Recovered', 'Total Deaths'};
legend(legendNames)
xlabel('Time (d)')
ylabel('# of people')
grid on

%% Plot totals by variant 
figure(2)

% Define compartment to plot from the reshapeInput output variable
plt_cmprtmnt = E_t + I_t;

% Sum all the infected per age group per variant
plt_cmprtmnt_var = squeeze(sum(reshape(plt_cmprtmnt,size(plt_cmprtmnt,1), GenP.nAgeGroups,[]),2));

% Plot the compartment defined pcer variant
hL = plot(simDates, plt_cmprtmnt_var, 'LineWidth', 1.5);
hold on
hP = plot(InputData.CasesAdj.Date, sumVar, 'o');

% Update color points
for j = 1:length(hP)
    hL(j).Color = colPoints(j,:);
    hP(j).Color = colPoints(j,:);
end

grid on
xlabel('Time (d)')
ylabel('# of people')
xlim([min(simDates), max(simDates)])
ylim([0 max(max(plt_cmprtmnt_var)*1.1)])
hold off
legend(EpiP.variantNames)

%% Plot percentage of circulating variant within population
figure(3)

% Plot the compartment defined per variant as percentage of circulating
% variants
plt_cmprtmnt_var_perc = plt_cmprtmnt_var ./ sum(plt_cmprtmnt_var,2) * 100; 
hA = area(simDates, plt_cmprtmnt_var_perc);
% hold on
% plot(InputData.CasesAdj.Collection_date, sumVar./sum(sumVar,2) * 100, 'o')
% hold off
grid on
xlabel('Time (d)')
ylabel('% of circulating variant')
ylim([0 100])
legend(EpiP.variantNames)
