function x_T = state_time(tabs, nro)

% t = tabs [T]
t = mod(tabs, nro.T);

%interpoled states
x=interp1(nro.tv,nro.yv(:,1),t,'spline');
y=interp1(nro.tv,nro.yv(:,2),t,'spline');
z=interp1(nro.tv,nro.yv(:,3),t,'spline');
xdot=interp1(nro.tv,nro.yv(:,4),t,'spline');
ydot=interp1(nro.tv,nro.yv(:,5),t,'spline');
zdot=interp1(nro.tv,nro.yv(:,6),t,'spline');

x_T=[x y z xdot ydot zdot];

end
