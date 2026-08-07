% This function calculates the constitutive component of the Cauchy
% stress tensor.

% theta = Fiber directions
% c = model coefficients
% F_ = deformation gradient
% N = Degree of polynomial model


function sigma_b = cauchy_stress_base(theta,c,F_,N)

a = c(1);
b = c(2);
c4 = c(3:end); % Coefficients for I4

S_ = [0,0,0;0,0,0;0,0,0];
C_ = transpose(F_)*F_;

for i=1:length(theta)
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    caux = c4(1:N-1); % Take only the coefficients for the fiber family
    dI4 = Inv4(C_,n,"der");
    S_ = S_ +...
            PolySEDFder_I4(n,caux,C_)*dI4;
    c4 = c4(N:end); % Update coefficients for next iteration
end
% We define the invariants I1 and I2
dI1 = Inv1(C_,"der");
dI2 = Inv2(C_,"der");

% We include the isotropic component
S_ = 2*(a*dI1 + b*dI2 + S_);
sigma_b = F_*S_*transpose(F_);

end
