% This function calculates the average relative error of the model 

% ODF = XRD data
% strain = strain data
% stress = stress data
% c = fitted coefficients of the model
% model = model used for the fitting
% N = Degree of polynomial for each anisotropic invariant
% M = Number of fiber families considered


function error = AvFitRelativeError(theta,strain,stress,c,model,N,M)

%n = length(strain);
ExpStress = EvalSQfit(theta,strain,c,model,N,M);

res = norm(ExpStress-stress);
aux = norm(stress);
error = (res/aux)*100;

end