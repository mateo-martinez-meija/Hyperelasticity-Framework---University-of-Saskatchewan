% This function calculates the average relative error of the model 

% strain = strain data
% stress = stress data
% c = fitted coefficients of the model
% model = model used for the fitting
% N = Degree of polynomial for each anisotropic invariant

function error = AvFitRelativeError(theta,strain,stress,c,model,N)

ExpStress = EvalSQfit(theta,strain,c,model,N);

res = norm(ExpStress-stress);
aux = norm(stress);
error = (res/aux)*100;

end
