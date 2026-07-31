% This function runs the commands to fit the polynomial model to the
% mechanical data, using the XRD data. It also generates the plot for fit.
% This function will create the plot for any SF sample.

% XRD_Data = Fiber orientation data from XRD analysis
% Mech_Data_Par = Mechanical data from stress applied parallel to fibers
% Mech_Data_Perp = Mechanical data from stres applied perpendicular to
%                  fibers
% n = Degree of polynomial for each anisotropic invariant
% m = Number of fiber families considered

function [coeff,Error] = CNC_SF_PolyFit(Mech_Data_Par,Mech_Data_Perp,n,m)

theta = FiberDirections(m);
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
fun_SF_Par = @(cPoly,xdata_Par)EvalStress3D_Full_I4(theta,xdata_Par,...
             cPoly,["MR","Poly"],n);
fun_SF_Perp = @(cPoly,xdata_Perp)EvalStress3D_Full_I4(theta,...
              xdata_Perp,cPoly,["MR_Perp","Poly"],n);
fun_SF = @(cPoly,strain)EvalSQfit(theta,strain,cPoly,["MR","Poly"],n);

% Initial coefficients
n = n-1; % Update n to be number of terms instead of degree of polynomial
aux1 = [10,10];
aux2 = ones(1,n*m);
cPoly0 = [aux1,aux2];

% Define the constraints for stability
eps = realmin;
aux3 = zeros(1,n*m);
A1 = [-1,0,aux3]; % First stability condition
A2 = [0,-1,aux3]; % Second stability condition
Aaux = [];
for i = 1:length(theta)
    aux4 = zeros(1,n*m);
    aux4(n*(i-1) + 1) = -1;
    Aaux = [Aaux;0,0,aux4]; % Additional stability conditions
end

A = [A1;A2;Aaux];
b = [-eps;-eps;-eps+zeros(m,1)];
lb = [];
ub = [];

% Least Square Curve Fit
coeff = lsqcurvefit(fun_SF,cPoly0,strain,stress,lb,ub,A,b);

%--------- ERROR -----------

Error = AvFitRelativeError(theta,strain,stress,coeff,["MR","Poly"],...
    n+1,m); % Change to n+1 since we updated n-->n-1 previously

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
    'Polynomial (' + string(n+1) +  ') model - ' + string(m) +...
    ' fiber families (SF \mid\mid)', ['SF distribution data' ...
    '(\perp)'],'Polynomial (' + string(n+1) +  ') model - ' + ...
    string(m) + ' fiber families (SF \perp)')
legend('Location','northwest')

end
