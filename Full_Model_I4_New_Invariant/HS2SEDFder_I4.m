% This function evaluates the derivative of the second Hargan and 
% Saccomandi strain energy density function for a given value of strain, 
% angle and set of coefficients.

% n = direction of the fiber
% c = coefficients
% C_ = Right Cauchy-Green tensor 

function eval = HS2SEDFder_I4(n,c,C_)

%C_ = [1/lambda,0,0;0,lambda^2,0;0,0,1/lambda];
%n = [cos(theta_i),sin(theta_i),0];
I4 = n*(C_*transpose(n));

eval = 0;
eval = eval + 2*c(1)*c(2)*(I4-1)/(c(2)-(I4-1)^2);

end