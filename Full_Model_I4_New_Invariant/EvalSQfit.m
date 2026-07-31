% This function evaluates the first half of the data using the EvalStress
% function for the strain parallel to the fibers and second half of the
% data using the EvalStress function for the strain perpendicular to the
% fibers. Currently the only option for the isotropic model is 'MR', which
% refers to the Mooney-Rivlin model.

% s = strain perentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF
% M_ = Structure tensor

function eval = EvalSQfit(M_,s,c,model)

iso = model(1);
aniso = model(2);
n = length(s);
xPar = s(1:n/2);
xPerp = s(1+n/2:end);

switch iso

    case "MR"
        ParStress = EvalStress3D_Full_I4(M_,xPar,c,["MR",aniso]);
        PerpStress = EvalStress3D_Full_I4(M_,xPerp,c,["MR_Perp",aniso]);

end
eval = [ParStress;PerpStress];
end
