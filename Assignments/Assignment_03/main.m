clear all; close all; clc;

calc_dynamics

% generate_trajectories
load('trajectories_mat.mat')


%% Question 2 - DESCENTRALIZED controller design

% Desired close loop performance:
zeta = 1;
wn = 20;

Bv =  subs(B, unkownvars, unkownvars_val);

Bbar = diag(diag(Bv));
for i=1:numel(Bbar)
    aux = children(vpa(Bbar(i)));
    Bbar(i) = sum(aux( ~has(aux,q) ));
end
Bbar = double(Bbar);

Km = diag(inv(fv));
Tm = diag(fv \ Bbar);

Tv = Tm;
Kv = 2*zeta*wn./Km;
Kp = wn^2./Km./Kv;


q0  = qs(1,2:4)';
dq0 = [0 0 0]';

try
    simdata = sim('omnibundle_control.slx','StopTime', num2str(qs(end,1)));
    simdata = simdata.simdata;
catch exception
    error(exception.message);
end


figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(simdata.time(:,1), simdata.signals(1).values(:,:), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'q measured'
legend({'q1','q2','q3'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(simdata.time(:,1), simdata.signals(2).values(:,:), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'q reference'
ax2.YLim = ax1.YLim;
ax3 = subplot(3,1,3);
plot(simdata.time(:,1), abs(simdata.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'error'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Decentralized configuration space control')
% fp.savefig('dec_control')


close all;
animate


%% Question 4 - CENTRALIZED controller design



q0  = qs(1,2:4)';
dq0 = [0 0 0]';
tau_d = [0;0;0];

try
    simdata = sim('omnibundle_control.slx','StopTime', num2str(qs(end,1)));
    simdata = simdata.simdata;
catch exception
    error(exception.message);
end


figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(simdata.time(:,1), simdata.signals(1).values, 'LineWidth',2 )
grid on, xlabel 't', ylabel 'q measured'
legend({'q1','q2','q3'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(simdata.time(:,1), simdata.signals(2).values, 'LineWidth',2 )
grid on, xlabel 't', ylabel 'q reference'
ax2.YLim = ax1.YLim;
ax3 = subplot(3,1,3);
plot(simdata.time(:,1), abs(simdata.signals(3).values), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'error'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Centralized inverse-dynamics control')
% fp.savefig('centr_control')

close all;
animate


%% Question 5 - Disturbance effect

syms d
tau_d = [0;d;0];
simplify(B \ tau_d)

q0  = qs(1,2:4)';
dq0 = [0 0 0]';
tau_d = [0;0.1;0];

try
    simdata = sim('omnibundle_control.slx','StopTime', num2str(qs(end,1)));
    simdata = simdata.simdata;
catch exception
    error(exception.message);
end


figure('Color','white','Position',[507   528   520   143]);
hold on;
plot(simdata.time(:,1), simdata.signals(1).values(:,1), 'LineWidth',2 )
plot(simdata.time(:,1), simdata.signals(2).values(:,1), 'LineWidth',2 )
plot(simdata.time(:,1), abs(simdata.signals(3).values(:,1)), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'q'
legend({'q measured','q reference','error'},'Location','northwest')
fp.savefig('centr_control_dist')




%% Question 7 -Operational Space Centralized Control

q0  = vpasolve(dir_kin_omni_bundle_fun(q)  - pes(1,2:4)');
q0 = double([q0.q1; q0.q2; q0.q3]);
dq0 = [0 0 0]';

tau_d = [0;0;0];


try
    aux = sim('omnibundle_control.slx','StopTime', num2str(pes(end,1)));
    simdata = aux.simdata;
    simdata_xe = aux.simdata_xe;
catch exception
    error(exception.message);
end


figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(simdata_xe.time(:,1), simdata_xe.signals(2).values(:,:), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'x_e measured'
legend({'x_e','y_e','z_e'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(simdata_xe.time(:,1), simdata_xe.signals(1).values(:,:), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'x_e reference'
ax2.YLim = ax1.YLim;
ax3 = subplot(3,1,3);
plot(simdata_xe.time(:,1), abs(simdata_xe.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't', ylabel 'error'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Centralized inverse-dynamics operational-space control')
% fp.savefig('centr_control_opsp')


qs = [0 0 0 0];
animate





