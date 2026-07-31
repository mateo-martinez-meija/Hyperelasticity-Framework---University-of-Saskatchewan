% This function evaluates the first half of the data using the EvalStress
% function for the strain parallel to the fibers and second half of the
% data using the EvalStress function for the strain perpendicular to the
% fibers. 

% ODF = XRD data
% s = strain perentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF
% N = Degree of polynomial for each anisotropic invariant
% M = Number of fiber families considered

function eval = EvalSQfit(theta,s,c,model,N,M)

iso = model(1);
aniso = model(2);
n = length(s);
xPar = s(1:n/2);
xPerp = s(1+n/2:end);

switch iso
    case "NH"
        ParStress = EvalStress3D_Full_I4(theta,xPar,c,["NH",aniso]);
        PerpStress = EvalStress3D_Full_I4(theta,xPerp,c,["NH_Perp",aniso]);

    case "MR"
        ParStress = EvalStress3D_Full_I4(theta,xPar,c,["MR",aniso],N,M);
        PerpStress = EvalStress3D_Full_I4(theta,xPerp,c,["MR_Perp",aniso],N,M);

end
eval = [ParStress;PerpStress];
end
