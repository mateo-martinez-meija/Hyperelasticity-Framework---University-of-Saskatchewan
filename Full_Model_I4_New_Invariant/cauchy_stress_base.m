% This function calculates the constitutive component of the Cauchy
% stress tensor.

% M_ = Structure tensor
% c = model coefficients
% F_ = deformation gradient

function sigma_b = cauchy_stress_base(M_,c,F_)

a = c(1);
b = c(2);
c4 = c(3:end); % Coefficients for I4

C_ = transpose(F_)*F_;

S_ = PolySEDFder_I4_star(C_,M_,c4)*M_;

% We define the invariants I1 and I2
dI1 = Inv1(C_,"der");
dI2 = Inv2(C_,"der");

% We include the isotropic component
S_ = 2*(a*dI1 + b*dI2 + S_);
sigma_b = F_*S_*transpose(F_);

end
