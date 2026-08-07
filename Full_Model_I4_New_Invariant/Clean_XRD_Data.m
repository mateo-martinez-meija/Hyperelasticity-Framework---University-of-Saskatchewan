% This function takes the original XRD data and returns a new dataset
% that is clean, symmetrized, and normalized. Additionally, the angles 
% are changed to radians. 

function new_data = Clean_XRD_Data(XRD_Data)

data = table2array(XRD_Data);
theta = data(:,1);
psi = data(:,2);

angle = transpose(linspace(2,178,1001));
eps = 0.75;

new_psi = [];
% Create a set of data with constant spacing between measurements
for i=1:numel(angle)
    lb = angle(i) - eps;
    ub = angle(i) + eps;

    pos = find(theta>lb & theta<ub);
    new_psi = [new_psi;mean(psi(pos))];
end

% Symmetrize data
sym_psi = zeros(numel(new_psi),1);
tol = 1e-6;
for j=1:(numel(angle)-1)/2
    pos = find(angle<(180 - angle(j)+tol) & angle>(180-angle(j)-tol));
    aux = (1/2)*(new_psi(j) + new_psi(pos));
    sym_psi(j) = aux;
    sym_psi(pos) = aux;
end
sym_psi(501) = new_psi(501);

% Change angle to radians and normalize
rad_angle = transpose(linspace(2*pi/180,178*pi/180,1001));
delta = rad_angle(2)-rad_angle(1);
Integral = delta*sum(sym_psi);
norm_sym_psi = sym_psi/Integral;

new_data = [rad_angle,norm_sym_psi];
%plot(rad_angle,norm_sym_psi,'b-')

end
