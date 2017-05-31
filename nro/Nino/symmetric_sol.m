function symm = symmetric_sol (Transfer)

% returns 1 if the angle between two vectors and the trajectory between its
% extreme points are in the same turn direction. 0 otherwise

aT = [0, Transfer(1,2) + Transfer(1,8), Transfer(1,3) + Transfer(1,9)];
bT = [0, Transfer(end,2) + Transfer(end,8), Transfer(end,3) + Transfer(end,9)];

P = cross(aT,bT);
S = sign(P(1)); 

if S < 0
   S = 0;
end

% if S == 1, angle is in the counterclockwise direction (angular velocity
% along x direction)

% ispolycw returns true if the polygonal contour vertices represented 
%by x and y are ordered in the clockwise direction. x and y are numeric vectors 
%with the same number of elements.

if S == ispolycw(Transfer(:,2)+ Transfer(:,8), Transfer(:,3) + Transfer(:,9)) 
    symm = 0;
else
    symm = 1;
end


