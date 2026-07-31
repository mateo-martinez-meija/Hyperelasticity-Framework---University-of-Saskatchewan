% We define a function for the second invariant

% C_ = right Cauchy-Green tensor
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv2(C_,type)

switch type

    case "inv"
        result = (1/2)*((trace(C_))^2 - trace(C_^2));

    case "der"
        result = trace(C_)*eye(length(C_)) - C_;

end
end
