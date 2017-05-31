A
x = sym('x');
y = sym('y');
z = sym('z');
xdot = sym('xdot');
ydot = sym('ydot');
zdot = sym('zdot');

x_T = sym('x_T');
y_T = sym('y_T');
z_T = sym('z_T');

mu = sym('mu');

% C1 = (1-mu)/d_T^3;
% C2 = mu/r_T^3;
% C3 = (1 - mu)/d^3;
% C4 = mu/r^3;
% d_T = sqrt( (x_T + mu)^2 + y_T^2 + z_T^2 );
% r_T = sqrt( (x_T - 1 + mu)^2 + y_T^2 + z_T^2 );
% d = sqrt( (x_T + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 );  %d_T + x
% r = sqrt( (x_T - 1 + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 ); %r_T + x

f = xdot;
g = ydot;
h = zdot;
l =   2*ydot + x + (1-mu)/sqrt( (x_T + mu)^2 + y_T^2 + z_T^2 )^3 * (x_T + mu) - (1 - mu)/sqrt( (x_T + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 )^3*(x_T + x + mu) + mu/sqrt( (x_T - 1 + mu)^2 + y_T^2 + z_T^2 )^3*(x_T + mu - 1) - mu/ sqrt( (x_T - 1 + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 )^3*(x_T + x + mu - 1);
m =  -2*xdot + y + (1-mu)/sqrt( (x_T + mu)^2 + y_T^2 + z_T^2 )^3 * y_T - (1 - mu)/sqrt( (x_T + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 )^3*(y + y_T) + mu/sqrt( (x_T - 1 + mu)^2 + y_T^2 + z_T^2 )^3*y_T - mu/ sqrt( (x_T - 1 + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 )^3*(y + y_T);
n =                (1-mu)/sqrt( (x_T + mu)^2 + y_T^2 + z_T^2 )^3 * z_T - (1 - mu)/sqrt( (x_T + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 )^3*(z + z_T) + mu/sqrt( (x_T - 1 + mu)^2 + y_T^2 + z_T^2 )^3*z_T - mu/ sqrt( (x_T - 1 + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 )^3*(z + z_T);

loc = [ f g h l m n];
var = [x y z xdot ydot zdot];

for i=1:6
    for j=1:6
        A(i,j) = diff(loc(i),var(j));
    end
end
