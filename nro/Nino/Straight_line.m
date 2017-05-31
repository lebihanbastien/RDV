function f = Straight_line(t,X)
f = zeros(6,1);


x_dot = X(4);
y_dot = X(5);
z_dot = X(6);


f(1) = x_dot;
f(2) = y_dot;
f(3) = z_dot;
%LVLH frame
f(4) = 0;
f(5) = 0;
f(6) = 0;




end