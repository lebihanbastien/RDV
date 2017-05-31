function f = Clohessy_Wiltshire(t,X,n)
f = zeros(6,1);


x = X(1);
y = X(2);
z = X(3);
x_dot = X(4);
y_dot = X(5);
z_dot = X(6);


f(1) = x_dot;
f(2) = y_dot;
f(3) = z_dot;
% %Orbital frame
% f(4) = 3 * n^2 * x + 2 * n * y_dot;
% f(5) = - 2 * n * x_dot;
% f(6) = - n^2 * z; 

%LVLH frame
f(4) = 2 * n * z_dot;
f(5) = - n^2 * y;
f(6) = 3 * n^2 * z - 2 * n * x_dot; 




end