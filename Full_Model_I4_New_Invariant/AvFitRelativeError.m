% This function calculates the average relative error of the model 

% ODF = XRD data
% strain = strain data
% stress = stress data
% c = fitted coefficients of the model
% model = model used for the fitting
% M_ = structure tensor


function error = AvFitRelativeError(M_,strain,stress,c, model)

%n = length(strain);
ExpStress = EvalSQfit(M_,strain,c,model);

res = norm(ExpStress-stress);
aux = norm(stress);
error = (res/aux)*100;

end