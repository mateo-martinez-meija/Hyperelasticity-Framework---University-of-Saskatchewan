% This function checks that the energy density function remains positive
% throughout the application of stress on the sample. Using the SEDFs
% defined previously, we will calculate the macro-energy of the system for
% each value of strain. 

function Energy = MREnergyCheck(ODF,s,c)

a = c(1);
b = c(2);
c4 = c(3:end); % Coefficients for I4
% p = c(end); % Hydrostatic pressure
lambda = 1 + s/100; % We calculate the strain ratio 
data = table2array(ODF);
theta = data(:,1); % XRD angle
psi = data(:,2); % XRD intensity

% We consider a stress applied parallel to the fibers which are pointing in
% the direction of the x-axis
F_ = [lambda,0,0;0,1/(lambda)^(1/2),0;0,0,1/lambda^(1/2)];
C_ = transpose(F_)*F_;

Energy = 0;

for i=1:length(psi)-1
    
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    delta = theta(i+1)-theta(i);

    Energy = Energy + PolySEDF(n,c4,C_)*psi(i)*delta;
end

I1 = trace(C_);
I2 = (1/2)*(trace(C_)^(2) - trace(C_^2));

Energy = Energy + a*(I1-3) + b*(I2-3);

end
