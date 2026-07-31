% This function calculates the fiber directions for the fiber families. We consider an odd number of families. If n
% is an even number then, we will use n+1 for the script. 

% n = number of fiber families

% theta = angles for fiber families

function theta = FiberDirections(n)

% Make n odd by adding 1 if it's even
if mod(n,2) == 0
    n = n+1;
end

theta = zeros(n,1);

theta(1) = pi/2; % Set first angle to pi/2
delta = pi/n; % Set the change in theta

k = (n-1)/2;

% Set the angles for the fiber families
for i=1:k
    theta(2*i) = pi/2 + i*delta;
    theta(2*i + 1) = pi/2 - i*delta;
end

end
