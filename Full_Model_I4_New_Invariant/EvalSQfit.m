% This function evaluates the first half of the data using the EvalStress
% function for the strain parallel to the fibers and second half of the
% data using the EvalStress function for the strain perpendicular to the
% fibers. Currently the only option for the isotropic model is 'MR', which
% refers to the Mooney-Rivlin model.

% s = strain perentages from mechanical data
% c = coefficients for the model
% M_ = Structure tensor

function eval = EvalSQfit(M_,s,c)

n = length(s);
xPar = s(1:n/2);
xPerp = s(1+n/2:end);

ParStress = EvalStress3D_Full_I4(M_,xPar,c,'parallel');
PerpStress = EvalStress3D_Full_I4(M_,xPerp,c,'perpendicular');

end
eval = [ParStress;PerpStress];
end
