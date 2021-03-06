function [zlci, alci] = syn2lci_acc(t, zsyn, asyn, mu)
% Change of coordinates: from Earth-Moon synodic to lunar-centered inertial
% coordinates. The acceleration vector is included
%--------------------------------------------------------------------------
% Moon state in EM synodical coordinates
%--------------------------------------------------------------------------
rsyn_m = [1-mu ; 0 ; 0 ];

%--------------------------------------------------------------------------
% Matrices
%--------------------------------------------------------------------------
ct  = cos(t);
st  = sin(t);
R   = [+ct -st 0; +st +ct 0 ; 0 0 1];
Rd  = [-st -ct 0; +ct -st 0 ; 0 0 0];
Rdd = [-ct +st 0; -st -ct 0 ; 0 0 0];

%--------------------------------------------------------------------------
% COC: state
%--------------------------------------------------------------------------
zlci = [R zeros(3) ; Rd R]*([zsyn(1:3) - rsyn_m ; zsyn(4:6)]);

%--------------------------------------------------------------------------
% COC: acceleration
%--------------------------------------------------------------------------
alci = Rdd * (zsyn(1:3) - rsyn_m) + 2*Rd*zsyn(4:6) + R*asyn;

end




