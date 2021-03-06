clear all; close all; clc;


calc_dynamics

% generate_trajectories
load('trajectories_mat.mat')


%% Question 2 - DESCENTRALIZED controller design

% Desired close loop performance:
zeta = 1;
wn = 20;
% wn = 8; % lab

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


tau_d = [0;0;0];
q0  = qs(1,2:4)';
dq0 = [0 0 0]';

PID_P = Kv + Kp.*Kv.*Tv;
PID_I = Kp.*Kv;
PID_D = Kv.*Tv;

try
    simdata = sim('omnibundle_control.slx','StopTime', num2str(qs(end,1)));
    simdata = simdata.simdata;
catch exception
    error(exception.message);
end

figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'joint angle q [rad]'
for i=1:3
    plot(simdata.time(:,1), simdata.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', ['q' num2str(i) ' desired'] )
    plot(simdata.time(:,1), simdata.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', ['q' num2str(i) ' measured'] )    
end
ylim([-1 4])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
plot(simdata.time(:,1), abs(simdata.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [rad]'
legend({'error q1','error q2','error q3'},'Location','northeastoutside')
sgtitle('Decentralized joint space control')
drawnow; ax2.Position(3) = ax1.Position(3);
% fp.savefig('dec_control')


% close all;
% animate

load('Lab_measurements/JS_decentralized_control.mat')

figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'joint angle q [rad]'
for i=1:3
    plot(labdata.time(:,1), labdata.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', ['q' num2str(i) ' desired'] )
    plot(labdata.time(:,1), labdata.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', ['q' num2str(i) ' measured'] )    
end
ylim([-1 4])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
plot(labdata.time(:,1), abs(labdata.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [rad]'
ylim([0 0.15])
legend({'error q1','error q2','error q3'},'Location','northeastoutside')
sgtitle('Decentralized joint space control - Laboratory')
drawnow; ax2.Position(3) = ax1.Position(3);
% fp.savefig('dec_control_lab')


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


figure('Color','white','Position',[96    54   479   346]); %[228   414   839   394]
subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'joint angle q [rad]'
for i=1:3
    plot(simdata.time(:,1), simdata.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', ['q' num2str(i) ' desired'] )
    plot(simdata.time(:,1), simdata.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', ['q' num2str(i) ' measured'] )    
end
ylim([-1 4])
% legend('Location','northeastoutside')
subplot(2,1,2);
plot(simdata.time(:,1), abs(simdata.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [rad]'
% legend({'error q1','error q2','error q3'},'Location','northeastoutside')
sgtitle('Centralized joint-space control')
% fp.savefig('centr_control')

% close all;
% animate



load('Lab_measurements/JS_centralized_control.mat')

figure('Color','white','Position', [604    45   593   346]);  %[228   414   839   394]
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'joint angle q [rad]'
for i=1:3
    plot(labdata.time(:,1), labdata.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', ['q' num2str(i) ' desired'] )
    plot(labdata.time(:,1), labdata.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', ['q' num2str(i) ' measured'] )    
end
ylim([-1 4])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
plot(labdata.time(:,1), abs(labdata.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [rad]'
ylim([0 0.10])
legend({'error q1','error q2','error q3'},'Location','northeastoutside')
sgtitle('Centralized joint space control - Laboratory')
drawnow; ax2.Position(3) = ax1.Position(3);
% fp.savefig('centr_control_lab')



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


figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'joint angle q [rad]'
for i=1:3
    plot(simdata.time(:,1), simdata.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', ['q' num2str(i) ' desired'] )
    plot(simdata.time(:,1), simdata.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', ['q' num2str(i) ' measured'] )    
end
ylim([-1 4])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
% hold on, grid on, xlabel 't', ylabel 'q measured'
plot(simdata.time(:,1), abs(simdata.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [rad]'
legend({'error q1','error q2','error q3'},'Location','northeastoutside')
sgtitle('Centralized joint-space control with load disturbance')
drawnow; ax2.Position(3) = ax1.Position(3);
% fp.savefig('centr_control_dist')



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


figure('Color','white','Position',[96    54   479   346]); % [228   414   839   394]
subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'end-effector position [m]'
for i=1:3
    plot(simdata_xe.time(:,1), simdata_xe.signals(1).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', [char('x'+i-1) ' desired'] )
    plot(simdata_xe.time(:,1), simdata_xe.signals(2).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', [char('x'+i-1) ' measured'] )    
end
ylim([-0.1 0.25])
% legend('Location','northeastoutside')
subplot(2,1,2);
plot(simdata_xe.time(:,1), abs(simdata_xe.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [m]'
% legend({'error x','error y','error z'},'Location','northeastoutside')
sgtitle('Centralized inverse-dynamics operational-space control')
% fp.savefig('centr_control_opsp')


% qs = [0 0 0 0];
% animate


load('Lab_measurements/OS_centralized_control_NF.mat')

figure('Color','white','Position',[604    45   593   346]); %[228   414   839   394]
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'end-effector position [m]'
for i=1:3
    plot(labdata_xe.time(:,1), labdata_xe.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', [char('x'+i-1) ' desired'] )
    plot(labdata_xe.time(:,1), labdata_xe.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', [char('x'+i-1) ' measured'] )    
end
ylim([-0.1 0.25])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
plot(labdata_xe.time(:,1), abs(labdata_xe.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [m]'
ylim([0 0.02])
legend({'error x','error y','error z'},'Location','northeastoutside')
drawnow; ax2.Position(3) = ax1.Position(3);
sgtitle('Centralized inverse-dynamics operational-space control')
% fp.savefig('centr_control_opsp_lab')





% % % % uiopen('Lab_experiments_data/OP_NF_Q2.fig',1)
% % % % dataObjs = findobj(gcf,'-property','DisplayName','Type','Line','-or','Type','Stair')
% % % % 
% % % % objq = findobj(dataObjs,'DisplayName','xe:1');
% % % % labdata_xe.signals(1).values(:,1) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe:2');
% % % % labdata_xe.signals(1).values(:,2) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe:3');
% % % % labdata_xe.signals(1).values(:,3) = objq.YData;
% % % % 
% % % % labdata_xe.time(:,1) = objq.XData;
% % % % 
% % % % objq = findobj(dataObjs,'DisplayName','xe_d:1');
% % % % labdata_xe.signals(2).values(:,1) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe_d:2');
% % % % labdata_xe.signals(2).values(:,2) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe_d:3');
% % % % labdata_xe.signals(2).values(:,3) = objq.YData;
% % % % 
% % % % objq = findobj(dataObjs,'DisplayName','e:1');
% % % % labdata_xe.signals(3).values(:,1) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','e:2');
% % % % labdata_xe.signals(3).values(:,2) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','e:3');
% % % % labdata_xe.signals(3).values(:,3) = objq.YData;
% % % % 
% % % % save('Lab_measurements/OS_centralized_control_NF','labdata_xe');





%% Question 9

syms hx

imp_force = simplify( [1 0 0] * J_geom(1:3,1:3) * inv(B) * J_geom(1:3,1:3).' * [hx 0 0].');

vpa(simplify(subs(imp_force, unkownvars, unkownvars_val)))



%% Question 10 -Operational Space Impedance Control

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


figure('Color','white','Position',[604    45   593   346]); %[228   414   839   394]
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'end-effector position [m]'
for i=1:3
    plot(simdata_xe.time(:,1), simdata_xe.signals(1).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', [char('x'+i-1) ' desired'] )
    plot(simdata_xe.time(:,1), simdata_xe.signals(2).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', [char('x'+i-1) ' measured'] )    
end
ylim([-0.1 0.25])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
plot(simdata_xe.time(:,1), abs(simdata_xe.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [m]'
legend({'error x','error y','error z'},'Location','northeastoutside')
sgtitle('Impedance operational-space control')
drawnow; ax2.Position(3) = ax1.Position(3);
% fp.savefig('impedance_cotntrol_opsp_OD')
% fp.savefig('impedance_cotntrol_opsp_UD')


% qs = [0 0 0 0];
% animate



%% Question 11


load('Lab_measurements/OS_impedance_UD.mat')

figure('Color','white','Position',[96    54   479   346]);
subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'end-effector position [m]'
for i=1:3
    plot(labdata_xe.time(:,1), labdata_xe.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', [char('x'+i-1) ' desired'] )
    plot(labdata_xe.time(:,1), labdata_xe.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', [char('x'+i-1) ' measured'] )    
end
ylim([-0.1 0.25])
% legend('Location','northeastoutside')
ax1 = subplot(2,1,2);
plot(labdata_xe.time(:,1), abs(labdata_xe.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [m]'
ylim([0 0.05])
% legend({'error x','error y','error z'},'Location','northeastoutside')
sgtitle('Impedance operational-space control')
linkaxes([ax1,ax2,ax3],'x')
% fp.savefig('impedance_cotntrol_opsp_UD_lab')


load('Lab_measurements/OS_impedance_OD.mat')

figure('Color','white','Position',[604    45   593   346]);
ax1 = subplot(2,1,1);
hold on, grid on, xlabel 't [s]', ylabel 'end-effector position [m]'
for i=1:3
    plot(labdata_xe.time(:,1), labdata_xe.signals(2).values(:,i), '-', 'LineWidth',3,'Color', fp.getColor(i,0.5),  'DisplayName', [char('x'+i-1) ' desired'] )
    plot(labdata_xe.time(:,1), labdata_xe.signals(1).values(:,i), '--', 'LineWidth',3,'Color', fp.getColor(i,1), 'DisplayName', [char('x'+i-1) ' measured'] )    
end
ylim([-0.1 0.25])
legend('Location','northeastoutside')
ax2 = subplot(2,1,2);
plot(labdata_xe.time(:,1), abs(labdata_xe.signals(3).values(:,:)), 'LineWidth',2 )
grid on, xlabel 't [s]', ylabel 'tracking error [m]'
ylim([0 0.05])
legend({'error x','error y','error z'},'Location','northeastoutside')
sgtitle('Impedance operational-space control')
drawnow; ax2.Position(3) = ax1.Position(3);
fp.savefig('impedance_cotntrol_opsp_OD_lab')


%%

% % % % uiopen('Lab_experiments_data/OP_F_T2.fig',1)
% % % % dataObjs = findobj(gcf,'-property','DisplayName','Type','Line','-or','Type','Stair')
% % % % 
% % % % objq = findobj(dataObjs,'DisplayName','xe:1');
% % % % labdata_xe.signals(1).values(:,1) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe:2');
% % % % labdata_xe.signals(1).values(:,2) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe:3');
% % % % labdata_xe.signals(1).values(:,3) = objq.YData;
% % % % 
% % % % labdata_xe.time(:,1) = objq.XData;
% % % % 
% % % % objq = findobj(dataObjs,'DisplayName','xe_d:1');
% % % % labdata_xe.signals(2).values(:,1) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe_d:2');
% % % % labdata_xe.signals(2).values(:,2) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','xe_d:3');
% % % % labdata_xe.signals(2).values(:,3) = objq.YData;
% % % % 
% % % % objq = findobj(dataObjs,'DisplayName','e:1');
% % % % labdata_xe.signals(3).values(:,1) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','e:2');
% % % % labdata_xe.signals(3).values(:,2) = objq.YData;
% % % % objq = findobj(dataObjs,'DisplayName','e:3');
% % % % labdata_xe.signals(3).values(:,3) = objq.YData;
% % % % 
% % % % save('Lab_measurements/OS_impedance_UD','labdata_xe');



