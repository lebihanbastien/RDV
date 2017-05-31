function fpa = flight_path_angle(yvc, center)
% FLIGHT_PATH_ANGLE gives the current flight path angle with respect to a 
% given center, in radians.
%
% If the state in position/velocity is (x, v),
% the sinus of the flight path angle with respect to the center x0
% is given by:
%
%   sin(fpa)  = - (x - xO).v / (|x - x0|*|v|)
%
% <a href="matlab: 
% web('https://en.wikipedia.org/wiki/Elliptic_orbit#Flight_path_angle','-browser')">(link)</a>
% 
%
% BLB 2016

% Position wrt object's center
erl = yvc(1:3) - center';

% Flight path angle is obtained from sinus
sinfpal = -  erl' * yvc(4:6)/ (norm(erl) * norm(yvc(4:6)));
fpa = asin(sinfpal);
end