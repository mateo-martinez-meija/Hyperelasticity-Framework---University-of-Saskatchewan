

function eigen = CheckPositiveDefinite(data)

data = table2array(data);
theta = data(:,1); % XRD angle
psi = data(:,2); % XRD intensity
M = [0,0,0;0,0,0;0,0,0];

for i=1:length(psi)-1
    
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    delta = theta(i+1)-theta(i);

    M = M + psi(i)*delta*transpose(n)*n;

[eigen_val,eigen_vec] = eig(M);

eigen = [eigen_val,eigen_vec];

end