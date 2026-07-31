% This function defines the custom SEDF that will allow us to calculate the
% total energy of the system. 

function Energy = CustSEDF(n,c,C_)

I4 = n*(C_*transpose(n));

Energy= (c(1)/(2*c(2)+1))*abs((I4-1)^(2*c(2)))*(I4-1);

end