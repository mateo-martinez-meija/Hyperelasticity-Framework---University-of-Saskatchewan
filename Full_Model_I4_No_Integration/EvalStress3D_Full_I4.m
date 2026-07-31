% Function EvalStress3D calculates the Cauchy Stress for each value of strain
% that we have in the Mechanical Data. The variable "model" lets us choose what
% type of Strain-Energy Density Function we want to use. Currently the only 
% model is the polynomic model ('Poly') with a Mooney-Rivlin isotropic 
% component.

% strain = array of strain percentages from mechanical data
% c = coefficients for the model
% model = array with model for isotropic and anisotropic SEDF
% N = Degree of polynomial for each anisotropic invariant

function Stress = EvalStress3D_Full_I4(theta,strain,c,model,N)

n = length(strain);
Iso = model(1);
Aniso = model(2);
Stress = [];
switch Iso
    case 'MR' % Mooney-Rivlin model
        for i=1:n
            Stress = [Stress;...
                %CauchyStress3D_Full_I4(data,strain(i),c4,Aniso)];
                MRCauchyStress3D_Full_I4(theta,strain(i),c,Aniso,N,M)];
        end
    case 'MR_Perp' % Mooney-Rivlin with perpendicular strain data
        for i=1:n
            Stress = [Stress;...
                MRCauchyStress3D_Full_I4_Perp(theta,strain(i),c,Aniso,N,M)];
        end
end
end
