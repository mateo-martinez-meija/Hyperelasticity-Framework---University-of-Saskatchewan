% This function calculates the structural tensor based on the Orientation
% Distribution Function

function M_ = StructTensor(psi,theta)

% Calculate the structural tensor
M_ = [0,0,0;0,0,0;0,0,0];
%M_ = [0,0;0,0];
for i=1:length(psi)-1
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    %n = [cos(theta(i)),sin(theta(i))];
    delta = theta(i+1)-theta(i);
    M_ = M_ + psi(i)*delta*transpose(n)*n;
end
end
