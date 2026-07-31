% This function evaluates the derivative of a custom strain energy 
% density function for a given value of strain, angle and set of 
% coefficients.

% n = direction of the fiber
% c = coefficients
% C_ = Right Cauchy-Green tensor

function eval = CustSEDFder_I4(n,c,C_)

%C_ = [1/lambda,0,0;0,lambda^2,0;0,0,1/lambda];
%n = [cos(theta_i),sin(theta_i),0];
I4 = n*(C_*transpose(n));

eval = c(1)*abs((I4-1)^(2*c(2)));

end