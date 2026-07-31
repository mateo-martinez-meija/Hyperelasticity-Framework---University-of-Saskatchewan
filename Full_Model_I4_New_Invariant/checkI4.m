% This function checks if I_4 becomes less than 1 at any point to see what
% should be the form of the logarithmic model.

function More1 = checkI4(ODF, strain)

m = length(strain);
ODF = table2array(ODF);
theta = ODF(:,1);
l = length(theta);
More1 = [];

for i=1:l
    n = [cos(theta(i)),sin(theta(i)),0]; % Fiber direction
    for j=1:m
        lambda = 1 + strain(j)/100; % We calculate the strain ratio
        F_ = [1/lambda^(1/2),0,0;0,lambda,0;0,0,1/lambda^(1/2)];
        C_ = transpose(F_)*F_;
        I4 = n*(C_*transpose(n));
        if I4>1
            More1 = [More1;[theta(i),strain(j),I4]];
        end
    end
end
