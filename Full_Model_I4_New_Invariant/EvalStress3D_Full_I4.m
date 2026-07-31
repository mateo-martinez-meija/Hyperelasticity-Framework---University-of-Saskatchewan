% Function EvalStress3D calculates the Cauchy Stress in the Lagrangian
% configuration for each value of strain that we have in the Mechanical
% Data. The variable "model" lets us choose what type of Strain-Energy
% Density Function we want to use. Currently the only models are the
% polynomic model ('Poly') and the exponential model ('Exp')

% data = XRD data of fiber orientation
% strain = array of strain percentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF
% N = Degree of polynomial for each anisotropic invariant
% M_ = Structure tensor

function Stress = EvalStress3D_Full_I4(M_,strain,c,model)

n = length(strain);
Iso = model(1);
Aniso = model(2);
Stress = [];
switch Iso
    case 'NH' % Neo-Hookean model
        %c4 = c(2:end);
        for i=1:n
            Stress = [Stress;...
                %CauchyStress3D_Full_I4(data,strain(i),c4,Aniso)];
                NHCauchyStress3D_Full_I4(data,strain(i),c,Aniso)];
        end
    case 'MR' % Mooney-Rivlin model
        %c4 = c(3:end);
        for i=1:n
            Stress = [Stress;...
                %CauchyStress3D_Full_I4(data,strain(i),c4,Aniso)];
                MRCauchyStress3D_Full_I4(M_,strain(i),c,Aniso)];
        end
    case 'NH_Perp' % Neo-Hookean with perpendicular strain data
        for i=1:n
            Stress = [Stress;...
                NHCauchyStress3D_Full_I4_Perp(data,strain(i),c,Aniso)];
        end
    case 'MR_Perp' % Mooney-Rivlin with perpendicular strain data
        for i=1:n
            Stress = [Stress;...
                MRCauchyStress3D_Full_I4_Perp(M_,strain(i),c,Aniso)];
        end
end
end