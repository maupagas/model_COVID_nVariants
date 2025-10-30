function [fImm] = calcImmunityLevel(fImm_0, nCases, nCasesMax, tpeak, tmaxdecay)

if nCases/nCasesMax > 1 || nCases/nCasesMax < 0
    fCases = 0;
else
    fCases = nCases/nCasesMax;
end

fImm = fImm_0 * (1-fCases) + (1-exp(-tpeak/tmaxdecay));