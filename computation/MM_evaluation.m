function MM = MM_evaluation(t_i, t_f, intervals, nro_T, options, cr3bp)

% divide TOF in sub-intervals
times = (t_f - t_i)/intervals;
t_int = t_i:times:t_f;

% variables
Phi_temp = zeros(6,6,length(t_int)-1);
phi0     = reshape(eye(6),1,[]);

for i=1:length(t_int)-1
%   target state for distances
    x0_T = reshape(state_time(t_int(i),nro_T),[],6);
   [~, abs_states] = ode113(@(t,y)single_abs(t,y,cr3bp.mu), [0 times], [x0_T phi0]', options); 
    Phi_temp(:,:,i) = reshape(abs_states(end,7:42),6,6);
end

% build the STM for the whole transfer
Phi = Phi_temp(:,:,end);

for i=1:length(t_int)-2
    Phi = Phi*Phi_temp(:,:,end-i);
end

MM = Phi;
