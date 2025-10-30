function [modelP, EpiP, GenP] = loadModelParam(filename)

% Load model parameters from each sheet 
modelP.nContacts    = readtable(filename, "Sheet", "nContacts");
modelP.pImm         = readtable(filename, "Sheet", "pImm");
modelP.infoVar      = readtable(filename, "Sheet", "infoVar");

% Epidemiological Parameters
EpiP_tbl = readtable(filename, "Sheet", "EpiP");
EpiP_tbl.Row = EpiP_tbl.Var1;
rowNames = EpiP_tbl.Var1;
EpiP_tbl.Var1 = [];

EpiP_M = table2array(EpiP_tbl);
% Assign values of each epidemiological parameter
for i = 1:length(rowNames)
    EpiP.(char(rowNames(i))) = EpiP_M(i,:);
end

% General Parameters
GenP_tbl = readtable(filename, "Sheet", "GenP");


for i = 1:length(GenP_tbl.Parameter)
    GenP.(char(GenP_tbl.Parameter(i))) = GenP_tbl.Values(i);
end
