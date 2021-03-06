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


%% Question 1

A_0_3_fun_aux = matlabFunction(A_0_3,'Vars',{q,L1,L2});
A_0_3_fun = @(q) A_0_3_fun_aux(q,0.132,0.132);

q_val = [0 0 0]';
A_0_3_fun(q_val)*[0 0 0 1]';

q_val = [0.67 -0.15 2.7]';
A_0_3_fun(q_val)*[0 0 0 1]';

q_val = [-0.73 0.25 1.5]';
A_0_3_fun(q_val)*[0 0 0 1]';


%% Question 2

% Match the position of the end-effector (T_0_3_fun(q)*[0 0 0 1]') in the
% base coordinates to a specific position p3

p3_val = [0.06 0.04 0.02 1]';
p3_resp = vpasolve( A_0_3_fun(q)*[0 0 0 1]' == p3_val );

fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(A_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]'))*[0 0 0 1]';

p3_val = [0.0847 0.0123 -0.005 1]';
p3_resp = vpasolve( A_0_3_fun(q)*[0 0 0 1]' == p3_val );

fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(A_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]'))*[0 0 0 1]';


%% Question 3

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


%% Question 2 (Assignment 02)

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

% fp.m2latex(pl{1})
% fp.m2latex(pl{2})
% fp.m2latex(pl{3})
% fp.m2latex(J_l{1})
% fp.m2latex(J_l{2})
% fp.m2latex(J_l{3})
% fp.m2latex(A_0_{1})
% fp.m2latex(A_0_{2})
% fp.m2latex(A_0_{3})
% fp.m2latex(T)
% fp.m2latex(U)


%% Question 3

% Equations of motion
gvec = jacobian(U,q).';

dB = sym(zeros(3,3));
for i=1:numel(B)
    dB(i) = jacobian(B(i),q)*dq;
end
coriolis_centr = dB*dq - 0.5*jacobian(dq.'*B*dq, q).';

for i=1:numel(q)
    H{i} = hessian(coriolis_centr(i),dq);
    simplify(dq.'*H{i}*dq - coriolis_centr(i))
end


f = diag(sym('f',[3 1]));
fs = diag(sym('fs',[3 1]));

DYNEQ = B*ddq + coriolis_centr + gvec - xi + f*dq + fs*sign(dq);

% gradL = gradient(L,dq);
% simplify(jacobian(gradL,q)*dq + jacobian(gradL,dq)*ddq - gradient(L,q) - (B*ddq + coriolis_centr + gvec))

% fp.m2latex(gvec)
% fp.m2latex(H{2}(2,3))
% fp.m2latex(H{3}(1,1))
% fp.m2latex(H{3}(1,2))
      
%% Question 4a

% simplify(DYNEQ)

% a = children(children(expand(DYNEQ(1))))
% subs(a{1}.', [q,dq,ddq], [[pi/6;pi/6;pi/6],[pi/6;pi/6;pi/6],[pi/6;pi/6;pi/6] ]) 

Pi = [I{1}(2,2); diag(I{2}); diag(I{3});
          L2^2*m3; L1*L2*m3; L2*m3;
          L1^2*m2; L1^2*m3; 
          L1*m3; L1*m2;
          diag(f);
          diag(fs)];
Pi_aux = sym('aux',[numel(Pi),1]);
DYNEQ_aux = subs(expand(DYNEQ), Pi, Pi_aux);

[Y,Tao]=equationsToMatrix(DYNEQ_aux,Pi_aux);
Y=vpa(simplify(Y));

if any(simplify( DYNEQ_aux - (Y*Pi_aux-Tao)))
    error('System identification failed!!')
end



% eqn2 = isolate(expand(DYNEQ(2))==0, L1*m2 )
% eqn2 = isolate(simplify(DYNEQ(2))==0, L1*(m3+0.5*m2) )
% solve(DYNEQ, L1*m3+0.5*L1*m2)

% Y(:,13), Y(:,14)
% Pi(13), Pi(14)
% Y(:,11), Y(:,12)
% Pi(11), Pi(12)

% Y(:,1), Y(:,10)
% Pi(1), Pi(7)

% for i=1:size(Y,2)
%     for j=i+1:size(Y,2)
%         fprintf('%f %f', i, j);
%         simplify(Y(:,i)./Y(:,j))
%         if simplify(Y(:,i)./Y(:,j)) == [1;1;1]
%             fprintf('%f %f', i, j);
%         end
%     end
% end

% 
% terms = children(expand(DYNEQ))
% terms{1}.'
if false
    for i=1:numel(Pi)
        fprintf('%s')
    end
end


if false
    for j=1:size(Y,1)
        for i=1:size(Y,2)
            if any(Y(:,i)~= 0)
                txt = strrep(latex(simplify(Y(j,i))), 'ddq', '\ddot q');
                txt = strrep(txt, 'dq', '\dot q');
                fprintf('\\scriptsize $%20s$ & \\scriptsize $%150s$ \\\\  \n',latex(Pi(i)),txt);
            end
        end
        fprintf('\n')
    end
    for j=1:size(Y,1)
        for i=1:size(Y,2)
            if any(Y(:,i)~= 0)
                fprintf('%20s * %s \n',Pi(i),simplify(Y(j,i)));
            end
        end
        fprintf('\n')
    end
end

%% Question 4a

resp = {};
for j=1:size(Y,2)
    if any(jacobian(Y(:,j),g) ~= 0)
        resp{1,end+1} = 0;
        for i=1:size(Y,1)
            resp{i,end} = Y(i,j);
        end
        fprintf('%s \\\\ \n',Pi(j))
    end
end


for i=1:size(resp,1)
    for j=1:size(resp,2)
        txt = strrep(latex(simplify(sym(resp{i,j}))), 'ddq', '\ddot q');
        txt = strrep(txt, 'dq', '\dot q');
        fprintf('%s ',txt);
        if j<size(resp,2)
            fprintf(' & ');
        end
    end
    fprintf(' & \\hdots \\\\ \n')
end


%% Question 4c
close all

% data = load('workspace_lab2_rev3.mat');

data = load('workspace_lab2_EXTRA_time.mat');
data = data.data;

figure('Color','white','Position',[257   380   854   437])
subplot(2,2,1)
plot(data.q(:,1), data.q(:,2:4)), grid on;
xlabel 'time [s]', ylabel 'Joint angles q [rad]';
legend({'q_1','q_2','q_3'})
subplot(2,2,2)
plot(data.Torque(:,1), data.Torque(:,2:4)), grid on;
xlabel 'time [s]', ylabel 'Joint torques \xi [Nm]';
legend({'\xi_1','\xi_2','\xi_3'})
ylim([-0.4 0.2]);
subplot(2,2,3)
plot(data.dq(:,1), data.dq(:,2:4)), grid on;
xlabel 'time [s]', ylabel 'Joint angles velocity  dq [rad/s]';
legend({'dq_1','dq_2','dq_3'})
ylim([-2 2]);
subplot(2,2,4)
plot(data.qdd(:,1), data.qdd(:,2:4)), grid on;
xlabel 'time [s]', ylabel 'Joint angles acceleration ddq [rad/s^2]';
legend({'ddq_1','ddq_2','ddq_3'})
ylim([-5 5]);
% fp.savefig('omni')


%% Question 4d


% % Data from the lab is not working. Intead, use simulation data from the TA
tmp = load('Data__Identification.mat');
data.q      = [tmp.sim_iden_q.time          tmp.sim_iden_q.signals.values];
data.dq     = [tmp.sim_iden_dq.time         tmp.sim_iden_dq.signals.values];
data.qdd    = [tmp.sim_iden_ddq.time        tmp.sim_iden_ddq.signals.values];
data.Torque = [tmp.sim_iden_torque.time     tmp.sim_iden_torque.signals.values];

% tmp = load('Data__Identification_2.mat');
% data.q      = [tmp.time  tmp.q];
% data.dq     = [tmp.time  tmp.dq];
% data.qdd    = [tmp.time  tmp.ddq];
% data.Torque = [tmp.time  tmp.Torque];



Y_fun = matlabFunction(Y,'Vars',{g q dq ddq});

% N = ceil(numel(data.time));
N = ceil(length(data.q));


n = size(Y,1);
Ydata = zeros(N*n, size(Y,2));
xidata = zeros(N*n, 1);
for i=1:N
    Ydata(n*i-n+1:n*i,:) = Y_fun(9.81, data.q(i,2:4)', data.dq(i,2:4)', data.qdd(i,2:4)');
    xidata(n*i-n+1:n*i,1) = data.Torque(i,2:4)';
end

Ybar = Ydata(1000:end,:);
xibar = xidata(1000:end,:);
% Pibar = (Ybar'*Ybar) \ Ybar' * xibar

idx = 1:20;
Pibar = pinv(Ybar(:,idx)'*Ybar(:,idx)) * Ybar(:,idx)' * xibar
% Pibar = Ybar(:,idx) \ xibar;

% Evaluate fitting performance
rms(Ybar*Pibar - xibar)



for i=1:numel(Pibar)
    fprintf('%20s: %8.4f \n',Pi(i),Pibar(i));
end



% Question 5

unkownvars = [I{1}(2,2); diag(I{2}); diag(I{3}); m2; m3; L1; L2; diag(f); diag(fs)];

idx = [1:7, 9, 11, 13,14, 15:20];
unkownvars_val_struct = vpasolve(Pibar(idx)==Pi(idx), unkownvars)
unkownvars_val = cellfun(@double, struct2cell(unkownvars_val_struct))

fnames = fieldnames(unkownvars_val_struct);
for i=1:numel(unkownvars_val)
    fprintf('%10s: %8.4f \n',fnames{i},unkownvars_val(i));
end


% unkownvars_val = [0.1  0 0.01 0.01  0 0.01 0.01  1 1  0.132 0.132  0.05 0.05 0.05 0.05 0.05 0.05]';

dyn_omni_bundle = subs(DYNEQ, unkownvars, unkownvars_val);
[Adyn,bdyn]=equationsToMatrix(dyn_omni_bundle,ddq);
dyn_omni_bundle_fun = matlabFunction( simplify(Adyn\bdyn), 'Vars', {g, q,dq,xi}, 'File', 'dyn_omni_bundle_fun');


% Simulate
torquein = data.Torque;
q0  = [0 0 0]';  % pi/2-0.1
dq0 = [1 0 0]';
Cfriction = 0;



try
    sim('omnibundle.slx')
catch
end



timeidx = simdata(:,1)<=10;
figure('Color','white','Position',[491   446   611   239])
plot(simdata(timeidx,1), simdata(timeidx,2:4)), grid on;
xlabel 'time [s]', ylabel 'Joint angles q [rad]';
legend({'q_1','q_2','q_3'})
% fp.savefig('sim')


% Animation

% simdata = data.q;
    
close all;

p1_fun = [0 0 0]';
p2_fun = matlabFunction( subs(p2(1:3),unkownvars, unkownvars_val), 'Vars', {q});
p3_fun = matlabFunction( subs(pe(1:3),unkownvars, unkownvars_val), 'Vars', {q});

figure('Color','white'), hold on, grid on, axis square;
lin2 = plot3( [0 [1 0 0]*p2_fun(q0)] ,[0 [0 1 0]*p2_fun(q0)], [0 [0 0 1]*p2_fun(q0)] ,'LineWidth',3);
lin3 = plot3( [0 [1 0 0]*p3_fun(q0)] ,[0 [0 1 0]*p3_fun(q0)], [0 [0 0 1]*p3_fun(q0)] ,'LineWidth',3);
xlim([-2*0.132 2*0.132])
ylim([-2*0.132 2*0.132])
zlim([-2*0.132 2*0.132])
view(45,15)
tic
for i=1:length(simdata)
    p2_val = p2_fun(simdata(i,2:4)');
    p3_val = p3_fun(simdata(i,2:4)');
    lin2.XData = [0 p2_val(1)];
    lin2.YData = [0 p2_val(2)];
    lin2.ZData = [0 p2_val(3)];
    lin3.XData = [p2_val(1) p3_val(1)];
    lin3.YData = [p2_val(2) p3_val(2)];
    lin3.ZData = [p2_val(3) p3_val(3)];
    drawnow();
    pause(simdata(i,1)-toc/2)
    if mod(i,10)==0
        fprintf('%f\n',toc)
    end
end

