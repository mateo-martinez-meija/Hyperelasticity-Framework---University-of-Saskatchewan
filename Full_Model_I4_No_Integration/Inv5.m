% We define a function for the fifth invariant

% C_ = right Cauchy-Green tensor
% n = normal vector with the direction of the fiber family
%       *Note: n is defined as a row vector
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv5(C_,n,type)

switch type

    case "inv"
        result = n*((C_^2)*transpose(n));

    case "der"
        result = transpose(n)*n*transpose(C_) + C_*(transpose(n))*n;

end
end
