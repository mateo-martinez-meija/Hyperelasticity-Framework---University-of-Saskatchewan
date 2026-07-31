% We define a function for the fourth invariant

% C_ = right Cauchy-Green tensor
% n = normal vector with the direction of the fiber family
%       *Note: n is defined as a row vector
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv4(C_,n,type)

switch type

    case "inv"
        result = n*(C_*transpose(n));

    case "der"
        result = transpose(n)*n;

end
end
