% This function evaluates the derivative of the mixed strain energy 
% density function for a given value of strain, angle and set of 
% coefficients. This function is made up of an quadratic anisotropic model
% for the SEDF up to a given value of strain and then a polynamial model
% for further values. 

% Note: The first entry of "c" has the coefficient for the quadratic model,
% while the other entries have the coefficients of the polynomial model.

% n = direction of the fiber
% c = coefficients
% C_ = Right Cauchy-Green tensor 
% s1 = strain value for where we stick both models together
% s2 = current value of strain where we evaluate the model

function eval = MixSEDFder_I4(n,s1,s2,c,C_)

%C_ = [1/lambda,0,0;0,lambda^2,0;0,0,1/lambda];
%n = [cos(theta_i),sin(theta_i),0];
I4 = Inv4(C_,n,"inv");

eval = 0;

if s2<s1
    eval = 2*c(1)*(I4-1);
else
    for i=1:length(c)-1
        j = i+1;
        eval = eval + (j)*c(i+1)*(I4-1)^(j-1);
    end
end

end