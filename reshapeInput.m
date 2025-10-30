function [S, E, I, R, D] = reshapeInput(Xin, nVar, nAgeGroups, inputFlag)
    
    % Determine the number of elements for each
    % compartment that is distributed by variant
    numElVar = nVar * nAgeGroups;

    % Identify each variable input first
    S = Xin(1:nAgeGroups, :);
    E = Xin(nAgeGroups + 1:(numElVar + nAgeGroups)*1,:);
    I = Xin(nAgeGroups + numElVar*1 + 1:nAgeGroups + numElVar*2, :);
    R = Xin(nAgeGroups + numElVar*2 + 1:nAgeGroups + numElVar*3, :);
    D = Xin(nAgeGroups + numElVar*3 + 1:nAgeGroups + numElVar*4, :);
    
%    % Check and correct for negative values
%     S(S < 0) = 0;   I(I < 0) = 0; 
%     R(R < 0) = 0;   D(D < 0) = 0; 
    % If this is the input for the model execution calculations, then
    % reshape. If it is the output over time, ignore
    if inputFlag
        % Reshape input to proper shape
        S = S';
        E = reshape(E, nAgeGroups, nVar)';
        I = reshape(I, nAgeGroups, nVar)';
        R = reshape(R, nAgeGroups, nVar)';
        D = reshape(D, nAgeGroups, nVar)';
    else
        S = S';
        E = E';
        I = I';
        R = R';
        D = D';
    end
end