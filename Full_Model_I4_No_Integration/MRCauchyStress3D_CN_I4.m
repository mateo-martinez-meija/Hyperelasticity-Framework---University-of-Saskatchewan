% This function calculates the Cauchy stress tensor in a 3 dimensional
% setup for a given value of strain. We start by calculating the second 
% Piola-Kirchhoff tensor numerically for a given strain energy density function
% for each individual fiber. We use the PK2 to calculate the Cauchy stress tensor.
% We consider a Mooney-Rivlin model for the isotropic component.

% s = strain percentage
% c = coefficients

function CauchyStress = MRCauchyStress3D_CN_I4(s,c)

a = c(1);
b = c(2);
c4 = c(4:end); % Coefficients for I4
lambda = 1 + s/100; % We calculate the strain ratio 
S_ = [0,0,0;0,0,0;0,0,0];

% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the y-axis
F_ = [1/lambda^(1/2),0,0;0,lambda,0;0,0,1/lambda^(1/2)];
C_ = transpose(F_)*F_;

% We define the derivatives for I1 and I2
dI1 = Inv1(C_,"der");
dI2 = Inv2(C_,"der");
I1 = Inv1(C_,"inv");
I2 = Inv2(C_,"inv");
%dI3 = Inv3(C_,"der");

% We include the isotropic component
S_ = 2*(a*dI1 + b*dI2);

% Find the Cauchy stress in the Lagrangian configuration
sigmaAux = (1/det(F_))*F_*S_*transpose(F_);
p = sigmaAux(3,3); % Calculate the hydrostatic pressure
sigma = sigmaAux - p*eye(3);

% We only consider sigma_yy to compare to experimental data
CauchyStress = sigma(2,2)/lambda; % We calculate the nominal stress
end
