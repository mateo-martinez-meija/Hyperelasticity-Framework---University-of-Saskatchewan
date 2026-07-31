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
% N = Degree of polynomial for each anisotropic invariant
% M = Number of fiber families considered

function CauchyStress = MRCauchyStress3D_Full_I4_Perp(theta,s,c,model,N,M)

a = c(1);
b = c(2);
%d = c(3);
c4 = c(3:end); % Coefficients for I4
% p = c(end); % Hydrostatic pressure
lambda = (1 + s/100)^2; % We calculate the strain ratio 
%data = table2array(data);
%theta = data(:,1); % XRD angle
%psi = data(:,2); % Weight for model
%[psi,theta] = weights(data,M);
S_ = [0,0,0;0,0,0;0,0,0];

c1 = [];
c2 = [];
m1 = [];
m2 = [];

for i=1:length(theta)
    c1 = [c1,c4(2*i-1)];
    c2 = [c2,c4(2*i)];
    m1 = [m1,(cos(theta(i)))^2];
    m2 = [m2,(sin(theta(i)))^2];
end

poly1 = 0;
poly2 = 0;
poly3 = a*lambda + b*lambda^2;
poly4 = 0;
poly5 =  - a - b*lambda;

for i=1:length(theta)
    poly1 = poly1 + 3*lambda*(c2(i)*m2(i)^3);
    poly2 = poly2 + 2*lambda*c1(i)*m2(i)^2 + ...
        6*lambda*c2(i)*(m2(i)^2)*(m1(i)*lambda - 1);
    poly3 = poly3 + 2*lambda*c1(i)*m2(i)*(m1(i)*lambda - 1) ...
        + 3*lambda*c2(i)*m2(i)*(m1(i)*lambda - 1)^2;
end

poly = [poly1,poly2,poly3,poly4,poly5];

aux = roots(poly);
%aux
%if sign(lambda1_fun(0))==sign(lambda1_fun(1))
%    CauchyStress = 0;
%    return
%else
%    lambda1 = fzero(lambda1_fun,x0)
j=4;
while ~isreal(aux(j)) || aux(j)<0
    j=j-1;
end
%if i~=4
%    aux
%end
lambda1 = aux(j);
lambda3 = 1/sqrt(lambda*lambda1);
lambda1 = sqrt(lambda1);
lambda = sqrt(lambda);


% We consider a stress applied perpendicular to the fibers which are pointing in
% the direction of the y-axis
% F_ = [lambda,0,0;0,1/(lambda)^(1/2),0;0,0,1/lambda^(1/2)];
F_ = diag([lambda,lambda1,lambda3]);
C_ = transpose(F_)*F_;

for i=1:length(theta)
    
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    caux = c4(1:N-1); % Taking the coefficients for the i-th invariant
    dI4 = Inv4(C_,n,"der"); % Define the derivative of I4

    % Numerical integration
    switch model
        case 'Poly'
            S_ = S_ +...
                PolySEDFder_I4(n,caux,C_)*dI4;
                %PolySEDFder_I4(n,caux,C_)*psi(i)*dI4;

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
    c4 = c4(N-1:end); % Update the coefficients for the next invariant
end

% We define the derivatives for I1 and I2
dI1 = Inv1(C_,"der");
dI2 = Inv2(C_,"der");
%dI3 = Inv3(C_,"der");

% We include the isotropic component
S_ = 2*(a*dI1 + b*dI2 + S_);

% Find the Cauchy stress in the Eulerian configuration
sigmaAux = (1/det(F_))*F_*S_*transpose(F_);
p = sigmaAux(3,3);
sigma = sigmaAux - p*eye(3);

% Using the push-forward mapping to find the stress in the Eulerian
% configuration
% aux = transpose(F_)\sigma/F_;

% We only consider sigma_xx to compare to experimental data
CauchyStress = sigma(1,1)/lambda;
end
