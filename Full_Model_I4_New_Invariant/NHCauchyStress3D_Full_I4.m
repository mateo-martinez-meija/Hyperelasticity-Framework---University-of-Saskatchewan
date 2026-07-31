% This function calculates the Cauchy stress tensor in a 3 dimensional
% Eulerian configuration setup for a given value of strain. We start by
% calculating the second Piola-Kirchhoff tensor numerically for a given
% strain energy density function for each individual fiber. We use the PK2
% to calculate the Cauchy stress tensor in the Lagrangian configuration and
% we push it forward by using the deformation gradient. We consider a
% Neo-Hookean model for the isotropic component.

% data = XRD data (angles and intensity)
% s = strain percentage
% c = coefficients
% model = SEDF for the fibers

function CauchyStress = NHCauchyStress3D_Full_I4(data,s,c,model)

mu = c(1);
c4 = c(2:end); % Coefficients for I4
% p = c(end); % Hydrostatic pressure
lambda = 1 + s/100; % We calculate the strain ratio 
data = table2array(data);
theta = data(:,1); % XRD angle
psi = data(:,2); % XRD intensity
S_ = [0,0,0;0,0,0;0,0,0];

% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the y-axis
F_ = [1/lambda^(1/2),0,0;0,lambda,0;0,0,1/lambda^(1/2)];
C_ = transpose(F_)*F_;

for i=1:length(psi)-1
    
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    delta = theta(i+1)-theta(i);

    % Numerical integration
    switch model
        case 'Poly'
            S_ = S_ +...
                PolySEDFder_I4(n,c4,C_)*psi(i)*delta*transpose(n)*n;
        
        case 'Exp'
            S_ = S_ +...
                ExpSEDFder_I4(n,c4,C_)*psi(i)*delta*transpose(n)*n;

        case 'HS1'
            S_ = S_ +...
                HS1SEDFder_I4(n,c4,C_)*psi(i)*delta*transpose(n)*n;

        case 'HS2'
            S_ = S_ +...
                HS2SEDFder_I4(n,c4,C_)*psi(i)*delta*transpose(n)*n;

        case 'Cust'
            S_ = S_ +...
                CustSEDFder_I4(n,c4,C_)*psi(i)*delta*transpose(n)*n;
    end
end

S_ = 2*((mu/2)*eye(3) + S_); % We include the isotropic component

% Find the Cauchy stress in the Lagrangian configuration
sigmaAux = (1/det(F_))*F_*S_*transpose(F_);
p = sigmaAux(3,3);
sigma = sigmaAux - p*eye(3);

% Using the push-forward mapping to find the stress in the Eulerian
% configuration
% aux = transpose(F_)\sigma/F_;

% We only consider sigma_yy to compare to experimental data
CauchyStress = sigma(2,2);
end