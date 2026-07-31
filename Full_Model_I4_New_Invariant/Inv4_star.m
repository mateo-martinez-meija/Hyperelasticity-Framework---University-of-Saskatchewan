% We define a function for the invariant I^{4*}

% C_ = right Cauchy-Green tensor
% M_ = Structure tensor
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv4_star(C_,M_,type)

switch type

    case "inv"
        result = trace(M_*C_);

    case "der"
        result = M_;

end
end
