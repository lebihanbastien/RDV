function [RC] = richardson_coefficients(cr3bp,li)
% RICHARDSON_COEFFICIENTS computation of the Richardson coefficients for
% third-order approximation of Halo orbits.
%
% RC = RICHARDSON_COEFFICIENTS(CR3BP, LI) stores in RC the Richardson
% coefficients for third-order approximation of Halo orbits around the
% Lagrange point LI of the system CR3BP. List of the coefficients:
%
% a21,a22, a23, a24, b21, b22, d21, a31, a32, b31, b32, d31, d32, s1, s2,
%  1   2    3    4    5    6    7    8    9    10   11   12   13  14   15
%
% l1, l2, lambda, omega_p, omega_v, kappa, d1, d2, c2, c3, c4,
% 16  17    18     19        20       21   22  23  24  25  26
%
% see: "Analytic construction of periodic orbits about the collinear points"
% <a href="matlab: 
% web('http://adsabs.harvard.edu/full/1980CeMec..22..241R','-browser')">(link)</a>
% Richardson 1980
%
% BLB 2015

%--------------------------------------------------------------------------
%Init
%--------------------------------------------------------------------------
mu = cr3bp.mu;

switch(li)
    
    case 1
        gamma_i = cr3bp.l1.gamma_i;
        
    case 2
        gamma_i = cr3bp.l2.gamma_i;
        
    case 3
        gamma_i = cr3bp.l3.gamma_i;
        
    otherwise
        gamma_i = cr3bp.l1.gamma_i;
        disp('gamma_i is initialized to its default value (l1-m2)');
end


%Cn coefficients for the corresponding li point
c2 = cn(mu, gamma_i, li, 2);
c3 = cn(mu, gamma_i, li, 3);
c4 = cn(mu, gamma_i, li, 4);

RC.c2 = c2;
RC.c3 = c3;
RC.c4 = c4;

%--------------------------------------------------------------------------
%Computation of characteristic constants of the periodic motion
%--------------------------------------------------------------------------
%lambda
lambda = sqrt(0.5*(2 - c2 + sqrt(9 * c2^2 - 8*c2)));
RC.lambda = lambda;
%omega_p
RC.omega_p = sqrt(0.5*(2 - c2 + sqrt(9 * c2^2 - 8*c2)));  %omega_p is set equal to lambda, see Koon et al. & Gomez et al. for details
%omega_v
RC.omega_v = sqrt(c2);
%kappa
kappa = (lambda^2+1+2*c2)/(2.0*lambda);
RC.kappa = kappa;
%Delta
RC.Delta = lambda^2-c2;

%Coeffs
d1 = (3*lambda^2)/kappa*( kappa*(6*lambda^2 - 1) - 2*lambda);
d2 = (8*lambda^2)/kappa*( kappa*(11*lambda^2 - 1) - 2*lambda);

RC.d1 = d1;
RC.d2 = d2;

%Calcul des coefficients
RC.d1 = (3*lambda^2)/kappa*( kappa*(6*lambda^2 - 1) - 2*lambda);

RC.d2 = (8*lambda^2)/kappa*( kappa*(11*lambda^2 - 1) - 2*lambda);

RC.a21 = (3*c3*(kappa^2 - 2))/(4*(1 + 2*c2));

RC.a22 = (3*c3)/(4*(1 + 2*c2));

RC.a23 = - (3*c3*lambda)/(4*kappa*d1)*(3*kappa^3*lambda - 6*kappa*(kappa-lambda) + 4);

RC.a24 = - (3*c3*lambda)/(4*kappa*d1)*(2 + 3*kappa*lambda);

RC.b21 = - (3*c3*lambda)/(2*d1)*(3*kappa*lambda - 4);

RC.b22 = (3*c3*lambda)/(d1);

RC.d21 = - c3/(2*lambda^2);

RC.a31 = - (9*lambda)/(4*d2)*(4*c3*(kappa*RC.a23 - RC.b21) + kappa*c4*(4 + kappa^2)) + (9*lambda^2 + 1 - c2)/(2*d2)*(3*c3*(2*RC.a23 - kappa*RC.b21) + c4*(2+3*kappa^2));

RC.a32 = - (9*lambda)/(4*d2)*(4*c3*(kappa*RC.a24 - RC.b22) + kappa*c4) - 3/(2*d2)*(9*lambda^2 + 1 - c2)*(c3*(kappa*RC.b22 + RC.d21 - 2*RC.a24) - c4);

RC.b31 = 3/(8*d2) * ( 8*lambda*(3*c3*(kappa*RC.b21 - 2*RC.a23) - c4*(2 + 3*kappa^2)) + (9*lambda^2 + 1 + 2*c2)*(4*c3*(kappa*RC.a23 - RC.b21) + kappa*c4*(4 + kappa^2)) );

RC.b32 =  1/d2    * ( 9*lambda*(c3*(kappa*RC.b22 + RC.d21 - 2*RC.a24) - c4) + 3/8*(9*lambda^2 + 1 + 2*c2)*(4*c3*(kappa*RC.a24 - RC.b22) + kappa*c4));

RC.d31 = 3/(64*lambda^2)*(4*c3*RC.a24 + c4);

RC.d32 = 3/(64*lambda^2)*(4*c3*(RC.a23 - RC.d21) + c4*(4 + kappa^2));

RC.s1 = (2*lambda*(lambda*(1 + kappa^2) - 2*kappa))^(-1 ) * (3/2*c3*(2*RC.a21*(kappa^2 - 2) -  RC.a23*(kappa^2 + 2) - 2*kappa*RC.b21)- 3/8*c4*(3*kappa^4 - 8*kappa^2 + 8));

RC.s2 = (2*lambda*(lambda*(1 + kappa^2) - 2*kappa))^(-1 ) * (3/2*c3*(2*RC.a22*(kappa^2 - 2) + RC.a24*(kappa^2 + 2) + 2*kappa*RC.b22 + 5*RC.d21) + 3/8*c4*(12 - kappa^2));

RC.l1 = - 3/2 * c3*(2*RC.a21 + RC.a23 +5*RC.d21) - 3/8*c4*(12 - kappa^2) + 2*lambda^2*RC.s1;

RC.l2 =   3/2 * c3*(RC.a24 - 2*RC.a22) + 9/8*c4 + 2*lambda^2*RC.s2;
end

%--------------------------------------------------------------------------
% cn parameters (see Richardson 1980)
%--------------------------------------------------------------------------
function cn = cn(mu, gamma_i, pointNumber, n)

switch(pointNumber)
    case 1
        cn =  gamma_i^(-3) * ( mu + (-1)^n*(1-mu)*gamma_i^(n+1)/(1-gamma_i)^(n+1) );
    case 2
        cn =  gamma_i^(-3) * ( (-1)^n*mu + (-1)^n*(1-mu)*gamma_i^(n+1)/(1+gamma_i)^(n+1) );
    case 3
        cn =  gamma_i^(-3) * ( 1 - mu + mu*gamma_i^(n+1)/(1+gamma_i)^(n+1) );
    otherwise
        cn = 0; %never here
end

end
