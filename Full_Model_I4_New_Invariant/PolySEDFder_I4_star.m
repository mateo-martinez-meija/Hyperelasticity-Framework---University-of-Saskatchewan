% This function evaluates the derivative of a polynomic strain energy 
% density function for a given value of strain, structure tensor and set of 
% coefficients.

% M_ = Structure tensor
% c = coefficients
% C_ = Right Cauchy-Green tensor 

function eval = PolySEDFder_I4_star(C_,M_,c)

I4 = Inv4_star(C_,M_,"inv");

eval = 0;
for i=1:length(c)
    j = i+1;
    eval = eval + (j)*c(i)*(I4-trace(M_))^(j-1);
end

end
