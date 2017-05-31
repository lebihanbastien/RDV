function [state,t] = theta2state(theta,nro)


t  = interp1(nro.alpha,nro.tv,theta,'spline');
state = state_time(t,nro);