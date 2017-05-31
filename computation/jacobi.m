function C = jacobi(yv, mu)
% JACOBI computes the jacobi constant.
%
% JACOBI(YV, MU) computes the jacobi constant at the state YV in the CRTBP
% system with the mass ratio MU.
%
% See equation 2.3.14 of Koon et al. 2006 <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
%
% BLB 2016

if(size(yv,2) == 3) %if only position is given, velocity is supposed null
    C = -2*u_barre(yv, mu);
else
    C = -2*u_barre(yv, mu) - (yv(4)^2 + yv(5)^2 +yv(6)^2);
end

end

function xout = u_barre(yv, mu)
%  Energy potential (from Koon et al. 2008)
%
mu1 = 1 - mu;
mu2 = mu;
r1 = sqrt( (yv(1)+mu2)^2 + yv(2)^2 + yv(3)^2 );
r2 = sqrt( (yv(1)-mu1)^2 + yv(2)^2 + yv(3)^2 );
xout = -1/2*(yv(1)^2+yv(2)^2) - mu1/r1 - mu2/r2 - 1/2*mu1*mu2;
end