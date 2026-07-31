% This function evaluates the derivative of a polynomic strain energy 
% density function for a given value of strain, fiber direction and set of 
% coefficients.

% n = direction of the fiber
% c = coefficients
% C_ = Right Cauchy-Green tensor 

function eval = PolySEDFder_I4(n,c,C_)

I4 = Inv4(C_,n,"inv");

eval = 0;
for i=1:length(c)
    j = i+1;
    eval = eval + (j)*c(i)*(I4-1)^(j-1);
end

end
