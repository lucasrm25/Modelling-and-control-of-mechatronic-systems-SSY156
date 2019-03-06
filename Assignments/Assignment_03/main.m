clear all; close all; clc;


%% OmniBundle parameters

% L1 = 0.132;
% L2 = 0.132;
syms L1 L2
q = sym('q',[3,1]);
dq = sym('dq',[3,1]);
ddq = sym('ddq',[3,1]);

DH_table = [0,-sym(pi)/2,0,q(1);
            L1,0,0,q(2);
            L2,0,0,q(3)-sym(pi)/2];
A = rb.T_DH_table(DH_table);        % A{i} = A^{i-1}_i

A_0_3 = A{1} * A{2} * A{3};


%% Calculate Kinematics

% 3 revolute joints
z0 = [0 0 1]';
z1 = A{1}(1:3,1:3)*z0;
z2 = A{1}(1:3,1:3)*A{2}(1:3,1:3)*z0;

p0 = [0 0 0 1]';            % selection of 4th column
p1 = A{1} * p0;
p2 = A{1} * A{2} * p0;
pe = A{1} * A{2} * A{3} * p0;

sel = [eye(3) zeros(3,1)];  % remove last row from vector

J_geom = [cross(z0,sel*(pe-p0)), cross(z1,sel*(pe-p1)), cross(z2,sel*(pe-p2));
          z0, z1, z2];
J_geom = simplify(J_geom);

syms phi theta psi
eulZYZ = [phi theta psi].';     % rotm2eul(T_0_3(1:3,1:3),'ZYZ')

J_analit = rb.T_zyz2geom(eulZYZ) \ J_geom;
J_analit = simplify(J_analit);


%% Calculate Dynamics

pl{1} = A{1} * p1;
pl{2} = A{1} * A{2} * [-L1/2 0 0 1].';
pl{3} = A{1} * A{2} * A{3} * [-L2/2 0 0 1].';

J_l{1} = [cross(z0,sel*(pl{1}-p0)), zeros(3,2);
          z0, zeros(3,2)];
J_l{2} = [cross(z0,sel*(pl{2}-p0)), cross(z1,sel*(pl{2}-p1)), zeros(3,1);
          z0, z1, zeros(3,1)];
J_l{3} = [cross(z0,sel*(pl{3}-p0)), cross(z1,sel*(pl{3}-p1)), cross(z2,sel*(pl{3}-p2));
          z0, z1, z2];

syms m2 m3
m = [0 m2 m3];
I{1} = diag(sym('I1_',[3,1]));
I{2} = diag(sym('I2_',[3,1]));
I{3} = diag(sym('I3_',[3,1]));


% A{i} = A^{i-1}_i
A_0_{1} = A{1};
A_0_{2} = A{1} * A{2};
A_0_{3} = A{1} * A{2} * A{3};

% Calculate Kinetic energy
B = zeros(3,3);
for i=1:3
    B = B + J_l{i}(1:3,:).' *                    m(i)                      * J_l{i}(1:3,:) + ...
            J_l{i}(4:6,:).' * A_0_{i}(1:3,1:3) * I{i} * A_0_{i}(1:3,1:3).' * J_l{i}(4:6,:);
end
B = simplify(B);
T = 0.5*dq.'*B*dq;
simplify(T);

% Calculate Potential energy
U = 0;
syms g
g0 = [0 0 -g].';    % gravity vector (g=9.8)
for i=1:3
    U = U - m(i) * g0.' *  pl{i}(1:3);
end
simplify(U);

% Generalized forces
xi = sym('xi',[3,1]);

L = T - U;


% Equations of motion
gvec = jacobian(U,q).';

dB = sym(zeros(3,3));
for i=1:numel(B)
    dB(i) = jacobian(B(i),q)*dq;
end
h = dB*dq - 0.5*jacobian(dq.'*B*dq, q).';

% for i=1:numel(q)
%     H{i} = hessian(h(i),dq);
%     simplify(dq.'*H{i}*dq/2 - h(i))
% end

C = jacobian(h,dq)/2;
% simplify(C*dq - h);

N = simplify(dB - 2*C);
% simplify(N+N.');          % Skew-symmetric

f = diag(sym('f',[3 1]));
fs = diag(sym('fs',[3 1]));

DYNEQ = B*ddq + h + gvec - xi + f*dq + fs*sign(dq);

gradL = gradient(L,dq);
simplify(jacobian(gradL,q)*dq + jacobian(gradL,dq)*ddq - gradient(L,q) - (B*ddq + h + gvec))

 

%Inertia body frame 1 
Iv{1} = [0.0036, 0.0031, 0.0045]; %Kg*m^2
Iv{2} = [0.0022, 0.0023, 0.0022]; %Kg*m^2
Iv{3} = [0.0010, 0.0011, 0.0009]; %Kg*m^2

%Masses
mv = [0.7, 0.11,0.08]; %Kg

%Lenghts 
L1v = 0.13; %m 
L2v = 0.13; %m

%Center of masses 
l1 = 0.07; %m 
l2 = 0.12; %m

%Friction (linear coefficients in the joint angular velocities) 
fv  = diag([0.0089, 0.0170, 0.0058]); %Nm*s/rad
fsv = diag([0, 0, 0]);

unkownvars     = [I{1}(2,2); diag(I{2}); diag(I{3}); m2; m3; L1; L2; diag(f); diag(fs)];
unkownvars_val = [Iv{1}(2); Iv{2}'; Iv{3}'; mv(2); mv(3); L1v; L2v; diag(fv); diag(fsv)];

dyn_omni_bundle = subs(DYNEQ, unkownvars, unkownvars_val);
[Adyn,bdyn]=equationsToMatrix(dyn_omni_bundle,ddq);
dyn_omni_bundle_fun = matlabFunction( simplify(Adyn\bdyn), 'Vars', {g, q,dq,xi}, 'File', 'dyn_omni_bundle_fun');

% test:
% dyn_omni_bundle_fun(9.8, q0, dq0, zeros(3,1))



%% Question 1 - Trajectory planning


% First two positions are ti and tf to next point. Next positions are joint angles
% qs = [ 0   0 0 0;
%        1   0 0 0;
%        3   0 pi/6 pi/6;
%        5   pi/4 0 pi/4;
%        7   0 pi/4 0;
%        9   pi/4 pi/4 pi/4];
   
% qs = [ 0  -0.6295 -0.2558 3.161      % 1   
%        2  -0.01275 -0.2342 2.993     % 5
%        3  -0.01275 -0.2342 2.993     % 5
%        5   0.361 -0.02856 2.364       % 9
%        6   0.361 -0.02856 2.364       % 9
%        8   0.6127 -0.2583 3.168       % 3
%        9   0.6127 -0.2583 3.168       % 3
%        11 -0.3619 -0.02611 2.363     % 7
%        12 -0.3619 -0.02611 2.363     % 7
%        14 -0.6295 -0.2558 3.161      % 1
%        15 -0.6295 -0.2558 3.161      % 1
%        ];
   
qs = [ 0  -0.6287 -0.2562 3.191      % 1
       2  -0.6287 -0.2562 3.191      % 1 
       4  -0.01028 -0.2411 3.023     % 5
       5  -0.01028 -0.2411 3.023     % 5
       7   0.3586 -0.03917 2.402     % 9
       8   0.3586 -0.03917 2.402     % 9
       10   0.602 -0.2652 3.197      % 3
       11   0.602 -0.2652 3.197      % 3
       13  -0.01028 -0.2411 3.023    % 5
       14  -0.01028 -0.2411 3.023    % 5
       16 -0.3614 -0.03509 2.393     % 7
       17 -0.3614 -0.03509 2.393     % 7
       19 -0.6287 -0.2562 3.191      % 1
       20 -0.6287 -0.2562 3.191      % 1
       ];
   
   
ddqc = 1;
dts = 0.05; 

[tv,qv,dqv,ddqv] = traj_planning(qs, ddqc,dts);


figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(tv,qv, 'LineWidth',2); 
grid on, xlabel 't', ylabel 'q'
legend({'q1','q2','q3'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(tv(1:end-1),dqv , 'LineWidth',2 );
grid on, xlabel 't', ylabel 'dq'
ax3 = subplot(3,1,3);
plot(tv(1:end-2), ddqv , 'LineWidth',2);
grid on, xlabel 't', ylabel 'ddq'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Configuration space trajectory planning')
% fp.savefig('trajectory_planning')



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


% fb_k = tf(zeros(numel(q)));
% for i=1:numel(q)
%     fb_k(i,i) = tf([Kv(i)*Tv(i) Kv(i)+Kp(i)*Kv(i)*Tv(i) Kp(i)*Kv(i)],[1 0])
% end


%% Simulate DESCENTRALIZED Control

q0  = qs(1,2:4)';
dq0 = [0 0 0]';

try
    simdata = sim('omnibundle_descentr_joint.slx','StopTime', num2str(qs(end,1)));
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



%% Question 4 - CENTRALIZED controller design

% dyn_omni_bundle = subs(DYNEQ, unkownvars, unkownvars_val);
[Adyn,bdyn]=equationsToMatrix(dyn_omni_bundle,xi);
inv_dyn_omni_bundle_fun = matlabFunction( simplify(Adyn\bdyn), 'Vars', {g, q,dq,ddq}, 'File', 'inv_dyn_omni_bundle_fun');


% test:
% a = inv_dyn_omni_bundle_fun(9.8, ones(3,1), ones(3,1), 10*ones(3,1))
% dyn_omni_bundle_fun(9.8, ones(3,1), ones(3,1), a)


zeta = 1;
wn = 20;

K_D = diag(repmat(2*zeta*wn,3,1));
K_P = diag(repmat(wn^2,3,1));



%% Simulate CENTRALIZED Control

q0  = qs(1,2:4)';
dq0 = [0 0 0]';
tau_d = [0;0;0];

try
    simdata = sim('omnibundle_centr_joint.slx','StopTime', num2str(qs(end,1)));
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
sgtitle('Centralized inverse-dynamics control')
% fp.savefig('centr_control')



%% Question 5 - Disturbance effect

syms d
tau_d = [0;d;0];
simplify(B \ tau_d)

q0  = qs(1,2:4)';
dq0 = [0 0 0]';
tau_d = [0;0.1;0];

try
    simdata = sim('omnibundle_centr_joint.slx','StopTime', num2str(qs(end,1)));
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


%% Operational Space Trajectory Planning


pes = [0   0.098 -0.078 -0.1;
       5   0.2 0 0;
       25   0.2 0 0];

ddpec = 0.5;
dts = 0.1;
[tv,pev,dpev,ddpev] = traj_planning(pes, ddpec, dts);

figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(tv,pev, 'LineWidth',2); 
grid on, xlabel 't', ylabel 'x_e'
legend({'x','y','z'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(tv(1:end-1),dpev , 'LineWidth',2 );
grid on, xlabel 't', ylabel 'dx_e'
ax3 = subplot(3,1,3);
plot(tv(1:end-2), ddpev , 'LineWidth',2);
grid on, xlabel 't', ylabel 'ddx_e'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Operational space trajectory planning')
fp.savefig('trajectory_planning_op_space')



%% Animate

% config
delay_factor = 2;       % 1 = no delay, 0.5 = half of real speed
fps = 60;

ta = 0:1/fps:simdata.time(end);
data_anim = [ta' interp1(simdata.time, simdata.signals(1).values, ta') ];


% figure('Color','white','Position',[491   446   611   239])
% plot(data_anim(:,1), data_anim(:,2:4)), grid on;
% xlabel 'time [s]', ylabel 'Joint angles q [rad]';
% legend({'q_1','q_2','q_3'})
% fp.savefig('sim')


% Animation

% data_anim = data.q;
    
% close all;
p1_fun = [0 0 0]';
p2_fun = matlabFunction( subs(p2(1:3),unkownvars, unkownvars_val), 'Vars', {q});
p3_fun = matlabFunction( subs(pe(1:3),unkownvars, unkownvars_val), 'Vars', {q});

figure('Color','white'), hold on, grid on, axis equal;

endeff_pos = p3_fun(qs(:,2:4)');
scatter3( endeff_pos(1,:), endeff_pos(2,:), endeff_pos(3,:), 'red' ,'filled');

lin2 = plot3( [0 [1 0 0]*p2_fun(q0)] ,[0 [0 1 0]*p2_fun(q0)], [0 [0 0 1]*p2_fun(q0)] ,'LineWidth',3);
lin3 = plot3( [0 [1 0 0]*p3_fun(q0)] ,[0 [0 1 0]*p3_fun(q0)], [0 [0 0 1]*p3_fun(q0)] ,'LineWidth',3);
xlim([0 2*0.132])
ylim([-0.15 0.15])
zlim([-0.15 0.1])
view(45,15)
tic
for i=2:length(data_anim)
    p2_val = p2_fun(data_anim(i,2:4)');
    p3_val = p3_fun(data_anim(i,2:4)');
    lin2.XData = [0 p2_val(1)];
    lin2.YData = [0 p2_val(2)];
    lin2.ZData = [0 p2_val(3)];
    lin3.XData = [p2_val(1) p3_val(1)];
    lin3.YData = [p2_val(2) p3_val(2)];
    lin3.ZData = [p2_val(3) p3_val(3)];
    
    plot3( [[1 0 0]*p3_fun(data_anim(i-1,2:4)') [1 0 0]*p3_fun(data_anim(i,2:4)')], ...
           [[0 1 0]*p3_fun(data_anim(i-1,2:4)') [0 1 0]*p3_fun(data_anim(i,2:4)')], ...
           [[0 0 1]*p3_fun(data_anim(i-1,2:4)') [0 0 1]*p3_fun(data_anim(i,2:4)')],...
           'LineWidth',3, 'Color', 'black');
    
    drawnow();
    pause(data_anim(i,1)-toc*delay_factor)
    if mod(i,10)==0
        fprintf('%f\n',data_anim(i,1))
    end
end


%% Help functions


function [tv,qv,dqv,ddqv] = traj_planning(qs, ddqc, dts)
    syms t
    qtcell = {};
    syms qt(t)
    for i=1:size(qs,1)-1
        for j=2:size(qs,2)
            qi = qs(i,j);
            qf = qs(i+1,j);
            dt = qs(i+1,1)-qs(i,1);
            if ddqc < 4*abs(qf-qi)/dt^2
                ddqc_n = 4*abs(qf-qi)/dt^2;
                fprintf('Max acceleration too low, changing to: %.2f\n',ddqc_n);
            else
                ddqc_n = ddqc;
            end
            ddqc_n = ddqc_n * sign(qf-qi);
            tc = dt/2 - 0.5*sqrt( (dt^2 *ddqc_n - 4*(qf-qi))/ddqc_n );
            if isnan(tc), tc=0; end

            qt(t) = piecewise( t<0,                  qi,...
                               (0<=t) & (t<=tc),     qi + 0.5*ddqc_n*t^2, ...
                               (tc<t) & (t<=dt-tc),  qi + ddqc_n*tc*(t-tc/2),...
                               (dt-tc<t) & (t<=dt),  qf - 0.5*ddqc_n*(dt-t)^2,...
                               t>dt,                 qf);
            qtcell{i,j-1} = @(t) double(qt(t-qs(i,1)));
        end
    end

    qv = [];
    tv=0:dts:qs(end,1);
    for i=1:numel(tv)
        idx =  min(find(tv(i) < qs(:,1)))-1;
        if isempty(idx)
            idx = size(qs,1)-1; 
        end
        qv(:,i) = cellfun(@(c) c(tv(i)),qtcell(idx,:),'UniformOutput',true)';
    end

    ddt = @(x,dt) conv2(x,[1 -1]/dt,'valid');

    dqv  =  ddt(qv,dts);
    ddqv = ddt(ddt(qv,dts),dts);
end















