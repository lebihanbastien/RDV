function [zlci] = syn2lci(t, zsyn, mu)
% Change of coordinates: from Earth-Moon synodic to lunar-centered inertial
% coordinates.
%--------------------------------------------------------------------------
% Moon state in EM synodical coordinates
%--------------------------------------------------------------------------
zsyn_m = [1-mu ; 0 ; 0 ;  0 ; 0 ; 0 ];

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
zlci = lciRotMat(t)*(zsyn - zsyn_m);

end


%--------------------------------------------------------------------------
% Subroutines
%--------------------------------------------------------------------------
function RotMat = lciRotMat(theta)
% Compute the rotation matrix associated to the angle theta
%
%   R =  | R11    0  |
%        | R21  R11  |
% with
%
%         | c -s 0 |          | -s -c 0 |
%   R11 = | s  c 0 |,   R21 = |  c -s 0 | 
%         | 0  0 1 |          |  0  0 0 |
% and
%       c = cos(theta), s = sin(theta)
%
% BLB 2016
RotMat11 = [+cos(theta) -sin(theta) 0; +sin(theta) +cos(theta) 0 ; 0 0 1];
RotMat21 = [-sin(theta) -cos(theta) 0; +cos(theta) -sin(theta) 0 ; 0 0 0];
RotMat   = [RotMat11 zeros(3) ; RotMat21 RotMat11];
end

