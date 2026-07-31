% Function EvalStress3D calculates the Cauchy Stress in the Lagrangian
% configuration for each value of strain that we have in the Mechanical
% Data. The variable "model" lets us choose what type of Strain-Energy
% Density Function we want to use. Currently the only models are the
% polynomic model ('Poly') and the exponential model ('Exp')

% data = XRD data of fiber orientation
% strain = array of strain percentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF
% s1 = strain value where the mixed model connects
% N = Degree of polynomial for each anisotropic invariant

function Stress = EvalStress3D_CN_I4(strain,c,model)

n = length(strain);
Stress = [];
switch model
    case 'NH' % Neo-Hookean model
        for i=1:n
            Stress = [Stress;...
                %CauchyStress3D_Full_I4(data,strain(i),c4,Aniso)];
                NHCauchyStress3D_Full_I4(data,strain(i),c)];
        end
    case 'MR' % Mooney-Rivlin model
        for i=1:n
            Stress = [Stress;...
                %CauchyStress3D_Full_I4(data,strain(i),c4,Aniso)];
                MRCauchyStress3D_CN_I4(strain(i),c)];
        end
end
end
