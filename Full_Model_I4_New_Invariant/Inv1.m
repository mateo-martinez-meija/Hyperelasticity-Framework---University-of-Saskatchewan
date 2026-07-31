% We define a function for the first invariant

% C_ = right Cauchy-Green tensor
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv1(C_,type)

switch type

    case "inv"
        result = trace(C_);

    case "der"
        result = eye(length(C_));

end
end
