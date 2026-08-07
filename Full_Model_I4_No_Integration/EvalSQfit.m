% This function evaluates the first half of the data using the EvalStress
% function for the strain parallel to the fibers and second half of the
% data using the EvalStress function for the strain perpendicular to the
% fibers. 

% theta = Fiber directions
% s = strain perentages from mechanical data
% c = coefficients for the model
% N = Degree of polynomial for each anisotropic invariant

function eval = EvalSQfit(theta,s,c,N)

n = length(s);
xPar = s(1:n/2);
xPerp = s(1+n/2:end);

ParStress = EvalStress3D_Full_I4(theta,xPar,c,'parallel',N);
PerpStress = EvalStress3D_Full_I4(theta,xPerp,c,'perpendicular',N);

eval = [ParStress;PerpStress];
end
