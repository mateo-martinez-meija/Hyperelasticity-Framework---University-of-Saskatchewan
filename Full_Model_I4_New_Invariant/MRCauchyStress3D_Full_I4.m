% This function calculates the Cauchy stress tensor in a 3 dimensional
% Eulerian configuration setup for a given value of strain. We start by
% calculating the second Piola-Kirchhoff tensor numerically for a given
% strain energy density function for each individual fiber. We use the PK2
% to calculate the Cauchy stress tensor in the Lagrangian configuration and
% we push it forward by using the deformation gradient. We consider a
% Mooney-Rivlin model for the isotropic component.

% data = XRD data (angles and intensity)
% s = strain percentage
% c = coefficients
% model = SEDF for the fibers
% M_ = Structure tensor

function CauchyStress = MRCauchyStress3D_Full_I4(M_,s,c,model)

a = c(1);
b = c(2);
c4 = c(3:end); % Coefficients for I4*
lambda = (1 + s/100)^2; % We calculate the stretch ratio squared
S_ = [0,0,0;0,0,0;0,0,0]; % Initialize the secon Piola-Kirchhoff tensor


c1 = c4(1);
c2 = c4(2);
m1 = M_(1,1);
m2 = M_(2,2);

% Calculate the stretches from the traction-free conditions
poly1 = 3*c2*(m1^3)*lambda;
poly2 = 2*c1*lambda*m1^2 + 6*c2*(m1^2)*lambda*(m2*(lambda-1)-m1);
poly3 = a*lambda + b*lambda^2 + 2*c1*m1*lambda*(m2*(lambda-1)-m1) ...
        + 3*c2*m1*lambda*(m2*(lambda-1)-m1)^2;
poly4 = 0;
poly5 = - a - b*lambda;
poly = [poly1,poly2,poly3,poly4,poly5];

aux = roots(poly);
% Choose the real positive root of the polynomial
i=4;
while ~isreal(aux(i)) || aux(i)<0
    i=i-1;
end
lambda1 = aux(i);

% Calculate the stretches
lambda3 = 1/sqrt(lambda*lambda1);
lambda1 = sqrt(lambda1);
lambda = sqrt(lambda);


% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the y-axis
F_ = diag([lambda1,lambda,lambda3]); % Deformation gradient
C_ = transpose(F_)*F_; % Right Cauchy-Green tensor

switch model
   case 'Poly'
            S_ = S_ +...
                PolySEDFder_I4_star(C_,M_,c4)*M_;
end

% We define the derivatives of I1 and I2
dI1 = Inv1(C_,"der");
dI2 = Inv2(C_,"der");


% We include the isotropic component
S_ = 2*(a*dI1 + b*dI2 + S_);

% Find the Cauchy stress
sigmaAux = (1/det(F_))*F_*S_*transpose(F_);
p = sigmaAux(3,3); % Calculate the hydrostatic pressure
sigma = sigmaAux - p*eye(3);


% We only consider sigma_yy to compare to experimental data
CauchyStress = sigma(2,2)/lambda; % We calculate the nominal stress
end
