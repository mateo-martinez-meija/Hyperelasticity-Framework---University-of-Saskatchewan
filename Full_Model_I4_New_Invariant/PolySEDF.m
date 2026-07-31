% This function defines the Polynomic SEDF that will allow us to calculate the
% total energy of the system. 

function Energy = PolySEDF(n,c,C_)

I4 = n*(C_*transpose(n));
Energy = 0;

for i=1:length(c)
    j = i+1;
    Energy = Energy + c(i)*(I4-1)^(j);
end

end