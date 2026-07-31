% This function calculates the weights for the each family of fibers based
% on the XRD normalized data. We consider an odd number of families. If n
% is an even number then, we will use n+1 for the script. 

% data = XRD data (angles and intensity)
% n = number of fiber families

% w = weights for new model
% theta = angles for fiber families

function [w,theta] = weights(data,n)

% Make n odd by adding 1 if it's even
if mod(n,2) == 0
    n = n+1;
end

if class(data) == "table"
    data = table2array(data);
end
angle = data(:,1); % XRD angle
psi = data(:,2); % Weight for model
eps = 0.001; % We define the length of the interval for average intensity

w = zeros(n,1);
theta = zeros(n,1);

w(1) = max(psi); % Set first weight for pi/2
theta(1) = pi/2; % Set first angle to pi/2
delta = pi/n; % Set the change in theta

k = (n-1)/2;

% Set the angles for the fiber families and respective weights
for i=1:k
    theta(2*i) = pi/2 + i*delta;
    theta(2*i + 1) = pi/2 - i*delta;
    
    if i==k
        break % We don't want to set w for 2k and 2k+1
    end

    % Set angles that define the needed interval for the average of psi
    % aux1 = (pi/2) - (i+1)*delta;
    % aux2 = (pi/2) - i*delta;
    % aux3 = (pi/2) + i*delta;
    % aux4 = (pi/2) + (i+1)*delta;

    %eps = 0.1; % We define the length of the interval
    ub1 = theta(2*i) + eps;
    ub2 = theta(2*i + 1) + eps;
    lb1 = theta(2*i) - eps;
    lb2 = theta(2*i+1) - eps;

    % Find the positions at which the intensities match the angles within
    % the intervals
    % pos = find((angle>aux1 & angle<aux2)|(angle>aux3 & angle<aux4));
    pos = find((angle>lb1 & angle<ub1)|(angle>lb2 & angle<ub2));

    temp = mean(psi(pos));
    w(2*i) = temp(end); % Weight is the average intensity
    w(2*i+1) = w(2*i); % Same weight due to symmetry
end

% Calculate the last weight to preserve normalization
aux5 = (1/2)*(1/delta - sum(w));
w(n-1) = aux5;
w(n) = aux5;

end