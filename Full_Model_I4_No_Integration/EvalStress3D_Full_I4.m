% Function EvalStress3D calculates the nominal Stress for each value of strain
% that we have in the Mechanical Data. Additionally, this function satisfies the
% the traction-free and shear-free conditions, and calculates the 
% appropriate stretches. 

% theta = Direction of each fiber family
% strain = array of strain percentages from mechanical data
% c = coefficients for the model
% loadcase = direction of load with respect to preferred fiber direction
% N = Degree of polynomial for each anisotropic invariant

function [Stress,residual] = EvalStress3D_Full_I4(theta,strain,c,loadcase,N)

n = length(strain);
u_prev = 0;
opts =  optimset('TolX',1e-12,'Display','off');

Stress = [];
switch loadcase
    case 'parallel' % Mooney-Rivlin model
        residual = [];
        for i=1:n
            lambda = 1 + strain(i)/100;
            fun = @(u)traction_residual(theta,exp(u),lambda,c,...
                loadcase,N);
            if abs(lambda-1)<1e-12
                u = 0;
            else
                u = fzero(fun,u_prev,opts);
            end
            lambda1 = exp(u);
            lambda3 = 1/(lambda*lambda1);
            F_ = diag([lambda1,lambda,lambda3]);

            sigma_b = cauchy_stress_base(theta,c,F_,N);
            p = sigma_b(3,3);
            sigma = sigma_b - p.*eye(3);

            P_axial = sigma(2,2)/lambda; % Calculate nominal stress
            Stress = [Stress;P_axial];
            residual = [residual;sigma(1,1),sigma(3,3),sigma(1,2)];

            u_prev = u;
        end
    case 'perpendicular' % Mooney-Rivlin with perpendicular strain data
        residual = [];
        for i=1:n
            lambda = 1 + strain(i)/100;
            fun = @(u)traction_residual(theta,exp(u),lambda,c,...
                loadcase,N);
            if abs(lambda-1)<1e-10
                u = 0;
            else
                u = fzero(fun,u_prev,opts);
            end
            lambda1 = exp(u);
            lambda3 = 1/(lambda*lambda1);
            F_ = diag([lambda,lambda1,lambda3]);

            sigma_b = cauchy_stress_base(theta,c,F_,N);
            p = sigma_b(3,3);
            sigma = sigma_b - p.*eye(3);

            P_axial = sigma(1,1)/lambda; % Calculate nominal stress
            Stress = [Stress;P_axial];
            residual = [residual;sigma(2,2),sigma(3,3),sigma(1,2)];

            u_prev = u;
        end
end



end
