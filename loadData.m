function [dataAD] = loadData(filename)

dataAD.Travelers = readtable(filename, "Sheet", "Travelers");
dataAD.TravelersM = table2array(dataAD.Travelers(:, 2:end));

dataAD.CasesAdj  = readtable(filename, "Sheet", "CasesAdj");


