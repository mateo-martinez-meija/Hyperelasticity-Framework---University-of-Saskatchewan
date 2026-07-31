% This function evaluates the integral over the angle theta for the
% anisotropic component to build the constraints on the new anisotropic 
% SEDF with the polynomial component to make it continuous up the second
% derivative. This integral is also evaluated for a specific value of
% strain "s", which is where we will 

% s = strain value at which we evaluate the integral
% k = exponent in the term (I4-1)^k
% ODF = Normalized XRD data for numerical integration
% orientation = fiber orientation, either 'par' or 'perp'

function eval = EvalIntPolyI4(s,k,ODF,orientation)

lambda = 1 + s/100; % We calculate the strain ratio
ODF = table2array(ODF);
theta = ODF(:,1); % XRD angle
psi = ODF(:,2); % XRD intensity

% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the y-axis

% par: stress applied parallel to the fibers which are pointing in the
% direction of the y-axis
% perp: stress applied perpendicular to the fibers which are pointing in
% the direction of the y-axis
switch orientation

    case "par"
        F_ = [1/lambda^(1/2),0,0;0,lambda,0;0,0,1/lambda^(1/2)];

    case "perp"
        F_ = [lambda,0,0;0,1/lambda^(1/2),0;0,0,1/lambda^(1/2)];

end
C_ = transpose(F_)*F_;
eval = 0;

for i=1:length(psi)-1

    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    delta = theta(i+1)-theta(i);

    % Numerical integration
    eval = eval + ((Inv4(C_,n,"inv")-1)^k)*psi(i)*delta;
end
end

