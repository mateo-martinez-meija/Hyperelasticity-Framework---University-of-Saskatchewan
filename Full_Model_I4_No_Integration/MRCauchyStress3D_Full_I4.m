% This function calculates the Cauchy stress tensor in a 3 dimensional
% setup for a given value of strain. We start by calculating the stretches 
% using the traction-free conditions. Then we calculate the second 
% Piola-Kirchhoff tensor numerically for a given
% strain energy density function. Currently, only using a polynomial SEDF. We 
% calcualte the hydrostatic pressure to find the Cauchy stress tensor. 
% Then, we calculate the nominal stress and return it. We consider a
% Mooney-Rivlin model for the isotropic component.

% theta = Fiber directions
% s = strain percentage
% c = coefficients
% model = SEDF for the fibers
% N = Degree of polynomial for each anisotropic invariant

function CauchyStress = MRCauchyStress3D_Full_I4(theta,s,c,model,N)

a = c(1);
b = c(2);
c4 = c(3:end); % Coefficients for I4 invariants
lambda = (1 + s/100)^2; % We calculate the stretch ratio squared 
S_ = [0,0,0;0,0,0;0,0,0]; % Initialize the second Piola-Kirchhoff tensor

c1 = [];
c2 = [];
m1 = [];
m2 = [];

% Define constants to help in the calculation of the stretches
for i=1:length(theta)
    c1 = [c1,c4(2*i-1)];
    c2 = [c2,c4(2*i)];
    m1 = [m1,(cos(theta(i)))^2];
    m2 = [m2,(sin(theta(i)))^2];
end

% Calculate the stretches from the traction-free conditions
poly1 = 0;
poly2 = 0;
poly3 = a*lambda + b*lambda^2;
poly4 = 0;
poly5 =  - a - b*lambda;

for i=1:length(theta)
    poly1 = poly1 + 3*lambda*(c2(i)*m1(i)^3);
    poly2 = poly2 + 2*lambda*c1(i)*m1(i)^2 + ...
        6*lambda*c2(i)*(m1(i)^2)*(m2(i)*lambda - 1);
    poly3 = poly3 + 2*lambda*c1(i)*m1(i)*(m2(i)*lambda - 1) ...
        + 3*lambda*c2(i)*m1(i)*(m2(i)*lambda - 1)^2;
end

% Define the polynomial for the traction-free boundary condition
poly = [poly1,poly2,poly3,poly4,poly5];

aux = roots(poly);

% Choose the real positive root of the polynomial
j=4;
while ~isreal(aux(j)) || aux(j)<0
    j=j-1;
end
%if i~=4
%    aux
%end
lambda1 = aux(j);

% Calculate the stretches
lambda3 = 1/sqrt(lambda*lambda1);
lambda1 = sqrt(lambda1);
lambda = sqrt(lambda);


% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the y-axis
F_ = diag([lambda1,lambda,lambda3]); % Deformation gradient
C_ = transpose(F_)*F_; % Right Cauchy-Green tensor

for i=1:length(theta)
    
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    caux = c4(1:N-1);
    dI4 = Inv4(C_,n,"der");

    % Numerical integration
    switch model
        case 'Poly'
            S_ = S_ +...
                PolySEDFder_I4(n,caux,C_)*dI4;
                %PolySEDFder_I4(n,caux,C_)*psi(i)*dI4;
    end

    c4 = c4(N-1:end); % Update the coefficients for the next invariant
end

% We define the derivatives for I1 and I2
dI1 = Inv1(C_,"der");
dI2 = Inv2(C_,"der");

% We include the isotropic component
S_ = 2*(a*dI1 + b*dI2 + S_);

% Find the Cauchy stress in the Lagrangian configuration
sigmaAux = (1/det(F_))*F_*S_*transpose(F_);
p = sigmaAux(3,3); % Calculate the hydrostatic pressure
sigma = sigmaAux - p*eye(3);

% We only consider sigma_yy to compare to experimental data
CauchyStress = sigma(2,2)/lambda; % We calculate the nominal stress
end
