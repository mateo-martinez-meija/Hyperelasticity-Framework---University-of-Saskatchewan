% This function calculates the average relative error of the model 

% strain = strain data
% stress = stress data
% c = fitted coefficients of the model
% M_ = structure tensor


function error = AvFitRelativeError(M_,strain,stress,c)

ExpStress = EvalSQfit(M_,strain,c);

res = norm(ExpStress-stress);
aux = norm(stress);
error = (res/aux)*100;

end
