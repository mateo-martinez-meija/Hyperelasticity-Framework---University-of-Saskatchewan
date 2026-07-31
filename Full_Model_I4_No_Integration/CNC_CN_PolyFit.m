% This function runs the commands to fit the polynomial model to the
% mechanical data, using the XRD data. It also generates the plot for fit.
% This function will create the plot for any CN sample.

% Mech_Data = Mechanical data

function [coeff,Error] = CNC_CN_PolyFit(Mech_Data)

Mech_Data = table2array(Mech_Data);

% Stress and strain data 
strain = Mech_Data(:,1);
stress = Mech_Data(:,2);

%------------ FIT ----------------

% Define the strain-energy density function 
fun_CN = @(c_MR,xdata_Par)EvalStress3D_CN_I4(xdata_Par,c_MR,"MR");

% Initial coefficients
cMR0 = [10,10];

% Define the constraints for stability
eps = realmin;
A1 = [-1,0]; % First stability condition
A2 = [0,-1]; % Second stability condition

A = [A1;A2]; 
b = [-eps;-eps];
lb = [];
ub = [];

% Least Square Curve Fit
coeff = lsqcurvefit(fun_CN,cMR0,strain,stress,lb,ub,A,b);

%--------- ERROR -----------

Error = AvFitRelativeErrorCN(strain,stress,coeff,"MR");

%--------- PLOT -----------

% Create a better set of x-values for the plot
x_full = linspace(0,max(strain),500);
% x_full = linspace(0,100,500);

% Plotting Curves
plot(strain,stress,'k.',x_full,fun_CN(coeff,x_full),'b-','MarkerSize', ...
    10,'LineWidth',5)
hold on
grid on
fontsize(16,"points")
xlabel('Strain $\varepsilon_{eng} \,\, (\%)$','Interpreter','latex')
ylabel("Stress (MPa)")
legend('CN distribution data','Mooney Rivlin model')
legend('Location','northwest')

end
