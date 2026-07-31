% We define a function for the third invariant

% C_ = right Cauchy-Green tensor
% type = etiher the invariant itself ("inv") or its derivative ("der")

function result = Inv3(C_,type)

switch type

    case "inv"
        result = det(C_);

    case "der"
        result = det(C_)*inv(C_);

end
end
