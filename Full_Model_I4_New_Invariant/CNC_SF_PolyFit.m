% This function runs the commands to fit the polynomial model to the
% mechanical data, using the XRD data. It also generates the plot for fit.
% This function will create the plot for any SF sample.

% XRD_Data = Fiber orientation data from XRD analysis
% Mech_Data_Par = Mechanical data from stress applied parallel to fibers
% Mech_Data_Perp = Mechanical data from stres applied perpendicular to
%                  fibers
% n = Degree of polynomial for each anisotropic invariant

function [coeff,Error] = CNC_SF_PolyFit(XRD_Data,Mech_Data_Par,Mech_Data_Perp,n)

data = table2array(XRD_Data);
theta = data(:,1); % XRD angle
psi = data(:,2); % XRD intensity
M_ = StructTensor(psi,theta);

Mech_Data_Par = table2array(Mech_Data_Par); % Parallel Stress
Mech_Data_Perp = table2array(Mech_Data_Perp); % Perpendicular Stress

% Stress and strain data from each dataset 
xdata_Par = Mech_Data_Par(:,1);
ydata_Par = Mech_Data_Par(:,2);
xdata_Perp = Mech_Data_Perp(:,1);
ydata_Perp = Mech_Data_Perp(:,2);

strain = [xdata_Par;xdata_Perp]; % Strain data for both tests
stress = [ydata_Par;ydata_Perp]; % Stress data for both tests

%------------ FIT ----------------

% Define the strain-energy density function 
fun_SF_Par = @(cPoly,xdata_Par)EvalStress3D_Full_I4(M_,xdata_Par,...
             cPoly,'parallel');
fun_SF_Perp = @(cPoly,xdata_Perp)EvalStress3D_Full_I4(M_,...
              xdata_Perp,cPoly,'perpendicular');
fun_SF = @(cPoly,strain)EvalSQfit(M_,strain,cPoly);

% Initial coefficients
n = n-1; % Update n to be number of terms instead of degree of polynomial
aux1 = [10,10];
aux2 = ones(1,n);
cPoly0 = [aux1,aux2];

% Define the constraints for stability
eps = realmin;
aux3 = zeros(1,n);
A1 = [-1,0,aux3]; % First stability condition
A2 = [0,-1,aux3]; % Second stability condition
aux3(1) = -1;
A3 = [0,0,aux3]; % Additional stability conditions

A = [A1;A2;A3];
b = [-eps;-eps;-eps];
lb = []; 
ub = [];

% Least Square Curve Fit
coeff = lsqcurvefit(fun_SF,cPoly0,strain,stress,lb,ub,A,b);

%--------- ERROR -----------

Error = AvFitRelativeError(M_,strain,stress,coeff);

%--------- PLOT -----------

% Create a better set of x-values for the plot
x_full = linspace(0,max(strain),500);

% Plotting Curves
plot(xdata_Par,ydata_Par,'b.',x_full,fun_SF_Par(coeff,x_full),'b-',...
    xdata_Perp,ydata_Perp,'r.',x_full,fun_SF_Perp(coeff,x_full),'r-',...
    'MarkerSize',10,'LineWidth',5)
hold on
grid on
fontsize(16,"points")
xlabel('Strain $\varepsilon_{eng} \,\, (\%)$','Interpreter','latex')
ylabel("Stress (MPa)")
legend('SF distribution data (\mid\mid)',...
    'Polynomial (' + string(n+1) + ...
    ') model -  ODF (SF \mid\mid)', ['SF distribution data' ...
    '(\perp)'],'Polynomial (' + string(n+1) +  ') model - ODF (SF \perp)')
legend('Location','northwest')

end
