% lambda1 = variable for traction-free condition
% lambda = stretch
% M_ = Structure tensor
% c = coefficients
% loadcase = direction of load (parallel or perpendicular)
% N = degree of polynomial for anisotropic component

function residual = traction_residual(M_,lambda1,lambda,c,loadcase)

isPar = strcmp(loadcase,'parallel');

lambda3 = 1/(lambda*lambda1);

if isPar
    F_ = diag([lambda1,lambda,lambda3]);
    sigma_b = cauchy_stress_base(M_,c,F_);

    residual = sigma_b(1,1) - sigma_b(3,3);

else
    F_ = diag([lambda,lambda1,lambda3]);
    sigma_b = cauchy_stress_base(M_,c,F_);

    residual = sigma_b(2,2) - sigma_b(3,3);
end


end
