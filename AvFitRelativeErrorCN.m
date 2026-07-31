% This function calculates the average relative error of the model 

% strain = strain data
% stress = stress data
% c = fitted coefficients of the model
% model = model used for the fitting

function error = AvFitRelativeErrorCN(strain,stress,c,model)

%n = length(strain);
ExpStress = EvalStress3D_CN_I4(strain,c,model);

res = norm(ExpStress-stress);
aux = norm(stress);
error = (res/aux)*100;

end