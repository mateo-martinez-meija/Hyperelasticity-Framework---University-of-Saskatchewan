% This function evaluates the derivative of a polynomic strain energy 
% density function for a given value of strain, angle and set of 
% coefficients.

% n = direction of the fiber
% c = coefficients
% C_ = Right Cauchy-Green tensor 

function eval = PolySEDFder_I4(n,c,C_)

%C_ = [1/lambda,0,0;0,lambda^2,0;0,0,1/lambda];
%n = [cos(theta_i),sin(theta_i),0];
I4 = Inv4(C_,n,"inv");

eval = 0;
for i=1:length(c)
    j = i+1;
    eval = eval + (j)*c(i)*(I4-1)^(j-1);
end

end
