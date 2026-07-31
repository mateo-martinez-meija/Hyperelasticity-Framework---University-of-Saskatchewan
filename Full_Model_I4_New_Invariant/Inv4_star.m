% We define a function for the fourth invariant

% C_ = right Cauchy-Green tensor
% n = normal vector with the direction of the fiber family
%       *Note: n is defined as a row vector
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv4_star(C_,M_,type)

switch type

    case "inv"
        result = trace(M_*C_);

    case "der"
        result = M;

end
end
