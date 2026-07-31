% Function EvalStress3D calculates the Cauchy Stress for each value of 
% strain that we have in the Mechanical Data. The variable "model" lets
% us choose what type of Strain-Energy Density Function we want to use.
% Currently the only model is the polynomic model ('Poly').

% strain = array of strain percentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF

function Stress = EvalStress3D_CN_I4(strain,c,model)

n = length(strain);
Stress = [];
switch model
    case 'MR' % Mooney-Rivlin model
        for i=1:n
            Stress = [Stress;...
                %CauchyStress3D_Full_I4(data,strain(i),c4,Aniso)];
                MRCauchyStress3D_CN_I4(strain(i),c)];
        end
end
end
