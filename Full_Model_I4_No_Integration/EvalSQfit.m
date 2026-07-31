% This function evaluates the first half of the data using the EvalStress
% function for the strain parallel to the fibers and second half of the
% data using the EvalStress function for the strain perpendicular to the
% fibers. 

% theta = Fiber directions
% s = strain perentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF
% N = Degree of polynomial for each anisotropic invariant

function eval = EvalSQfit(theta,s,c,model,N)

iso = model(1);
aniso = model(2);
n = length(s);
xPar = s(1:n/2);
xPerp = s(1+n/2:end);

switch iso
    case "MR"
        ParStress = EvalStress3D_Full_I4(theta,xPar,c,["MR",aniso],N);
        PerpStress = EvalStress3D_Full_I4(theta,xPerp,c,["MR_Perp",aniso],N);
end
eval = [ParStress;PerpStress];
end
