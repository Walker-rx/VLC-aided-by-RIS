clear 
close all
syms x y;

fh = @(x, y, z) x + y + z;
[numVars, varNames] = getFunctionHandleParams(fh);

disp(['Number of variables: ', num2str(numVars)]);
disp(['Variable names: ', varNames{:}]);