% This function check that the total energy is positive for every value of
% strain in the mechanical data. 

function Result = EvalEnergyCheck(ODF,strain,c)

n = length(strain);
% Iso = model(1);
% Aniso = model(2);
Energy = [];
Result = "The model checks out! :)";

for i=1:n
    NewEnergy = MREnergyCheck(ODF,strain(i),c);
    Energy = [Energy;NewEnergy];
    if NewEnergy<-1e-20
        Result = "The energy is not positive! T_T";
        %break
    end

end

end