% This function calculates the average relative error for any function
% using the fitted coefficients
% M_ = structure tensor

function error = AvRelativeError(M_,strain,stress,c, model)

%n = length(strain);
ExpStress = EvalStress3D_Full_I4(M_,strain,c,model);

res = norm(ExpStress-stress);
aux = norm(stress);
error = (res/aux)*100;

end