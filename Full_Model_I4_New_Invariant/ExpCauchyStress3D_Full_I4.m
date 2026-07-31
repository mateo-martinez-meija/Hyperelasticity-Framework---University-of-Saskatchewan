% This function calculates the Cauchy stress tensor in a 3 dimensional
% Eulerian configuration setup for a given value of strain. We start by
% calculating the second Piola-Kirchhoff tensor numerically for an
% exponential strain energy density function for each individual fiber. We 
% use the PK2 to calculate the Cauchy stress tensor in the Lagrangian 
% configuration and we push it forward by using the deformation gradient. 

% data = XRD data (angles and intensity)
% s = strain percentage
% c = coefficients

function CauchyStress = ExpCauchyStress3D_Full_I4(data,s,c)

mu = c(1);
c4 = c(2:end); % Coefficients for I4
lambda = 1 + s/100; % We calculate the strain ratio 
data = table2array(data);
theta = data(:,1); % XRD angle
psi = data(:,2); % XRD intensity
S_ = [0,0,0;0,0,0;0,0,0];

% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the y-axis
F_ = [1/lambda^(1/2),0,0;0,lambda,0;0,0,1/lambda^(1/2)];

for i=1:length(psi)-1
    
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    delta = theta(i+1)-theta(i);

    % Numerical integration
    S_ = S_ + ExpSEDFder_I4(lambda,n,c4)*psi(i)*delta*transpose(n)*n;

end

S_ = 2*((mu/2)*eye(3) + S_); % We include the isotropic component

% Find the Cauchy stress in the Lagrangian configuration
sigma = (1/det(det(F_)))*F_*S_*transpose(F_);

% Using the push-forward mapping to find the stress in the Eulerian
% configuration
aux = transpose(F_)\sigma/F_;

% We only consider sigma_yy to compare to experimental data
CauchyStress = aux(2,2);
end