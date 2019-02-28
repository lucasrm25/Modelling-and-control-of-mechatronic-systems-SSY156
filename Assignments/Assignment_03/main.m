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
coriolis_centr = dB*dq - 0.5*jacobian(dq.'*B*dq, q).';

for i=1:numel(q)
    H{i} = hessian(coriolis_centr(i),dq);
%     simplify(dq.'*H{i}*dq/2 - coriolis_centr(i));
end


f = diag(sym('f',[3 1]));
fs = diag(sym('fs',[3 1]));

DYNEQ = B*ddq + coriolis_centr + gvec - xi + f*dq + fs*sign(dq);

gradL = gradient(L,dq);
simplify(jacobian(gradL,q)*dq + jacobian(gradL,dq)*ddq - gradient(L,q) - (B*ddq + coriolis_centr + gvec))

 

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
fv  = [0.0089, 0.0170, 0.0058]; %Nm*s/rad
fsv = [0, 0, 0];

unkownvars     = [I{1}(2,2); diag(I{2}); diag(I{3}); m2; m3; L1; L2; diag(f); diag(fs)];
unkownvars_val = [Iv{1}(2); Iv{2}'; Iv{3}'; mv(2); mv(3); L1v; L2v; fv'; fsv'];

dyn_omni_bundle = subs(DYNEQ, unkownvars, unkownvars_val);
[Adyn,bdyn]=equationsToMatrix(dyn_omni_bundle,ddq);
dyn_omni_bundle_fun = matlabFunction( simplify(Adyn\bdyn), 'Vars', {g, q,dq,xi}, 'File', 'dyn_omni_bundle_fun');

% ddq = dyn_omni_bundle_fun(9.8, q0, dq0, dq0)


%% Simulate
torquein = [0 0 0 0; 1e10 0 0 0];
q0  = [0 0 0]';  % pi/2-0.1
dq0 = [1 0 0]';

try
    sim('omnibundle.slx')
catch exception
    error(exception.message);
end

%% Animate

timeidx = simdata(:,1)<=10;
figure('Color','white','Position',[491   446   611   239])
plot(simdata(timeidx,1), simdata(timeidx,2:4)), grid on;
xlabel 'time [s]', ylabel 'Joint angles q [rad]';
legend({'q_1','q_2','q_3'})
% fp.savefig('sim')


% Animation

% simdata = data.q;
    
% close all;
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



%% Question 1


% First two positions are ti and tf to next point. Next positions are joint angles
qs = [ 0  0 0 0;
       3  0 pi/6 pi/6;
       3.5  0 pi/6 pi/6;
       6  pi/4 0 pi/4;
       6.5  pi/4 0 pi/4;
       9  0 pi/4 0;
       9.5  0 pi/4 0;
       13  pi/4 pi/4 pi/4
       16  pi/4 pi/4 pi/4];

% Define trajectory equation: q(t) = a5*t^5 + ... + a1*t + a0
% a = sym('a',[3,6]);
% syms t
% qt  = a(:,6)*t^5 + a(:,5)*t^4 + a(:,4)*t^3 + a(:,3)*t^2 + a(:,2)*t^1 + a(:,1)*t^0;
% dqt = diff(qt,t,1);
% ddqt = diff(qt,t,2);
% params = solve( [  subs(qt,  t,trg(i,  1)) == trg(i,3:end)';...
%                    subs(qt,  t,trg(i+1,1))  == trg(i+1,3:end)';...
%                    subs(dqt, t,trg(i,  1))  == zeros(3,1);...
%                    subs(dqt, t,trg(i+1,1))  == zeros(3,1);...
%                    subs(ddqt,t,trg(i,  1))  == zeros(3,1);...
%                    subs(ddqt,t,trg(i+1,1))  == zeros(3,1)],a(:));
% params = cellfun(@double, struct2cell(params));

ddqc = 5;

syms t
qtcell = {};
for i=1:size(qs,1)-1
    for j=2:size(qs,2)
        qi = qs(i,j);
        qf = qs(i+1,j);
        dt = qs(i+1,1)-qs(i,1);
        tc = dt/2 - 0.5*sqrt( (dt^2 *ddqc - 4*(qf-qi))/ddqc );
               
        if ddqc < 4*abs(qf-qi)/dt^2
            ddqc_n = 4*abs(qf-qi)/dt^2;
            printf('Max acceleration too low, changing to: %.2f',ddqc);
        else
            ddqc_n = ddqc;
        end
        
        qt(t) = piecewise( t<=0,        0,...
                           0<t<=tc,     qi + 0.5*ddqc_n*t^2, ...
                           tc<t<=dt-tc, qi + ddqc_n*tc*(t-tc/2),...
                           dt-tc<t<=dt, qf - 0.5*ddqc_n*(dt-t)^2,...
                           dt<t,        0);
        qtcell{i,j-1} = @(t) double(qt(t-qs(i,1)));
    end
end


dt = 0.01;

qv = [];
t=dt:dt:5;
for i=1:numel(t)
    qv(:,i) = sum(cellfun(@(c) c(t(i)),qtcell,'UniformOutput',true))';
end

ddt = @(x,dt) conv2(x,[1 -1]/dt,'valid');

(qv(:,end) - qv(:,end-1))/dt

figure('Color','white');
hold on, grid on;
subplot(3,1,1)
plot(t,qv, 'LineWidth',2)
subplot(3,1,2)
plot(t(1:end-1), ddt(qv,dt) , 'LineWidth',2 )
subplot(3,1,3)
plot(t(1:end-2), ddt(ddt(qv,dt),dt) , 'LineWidth',2)


% A_0_3_fun_aux = matlabFunction(A_0_3,'Vars',{q,L1,L2});
% A_0_3_fun = @(q) A_0_3_fun_aux(q,0.132,0.132);
% 
% q_val = [0 0 0]';
% A_0_3_fun(q_val)*[0 0 0 1]';
% 
% q_val = [0.67 -0.15 2.7]';
% A_0_3_fun(q_val)*[0 0 0 1]';
% 
% q_val = [-0.73 0.25 1.5]';
% A_0_3_fun(q_val)*[0 0 0 1]';


%% Question 2

% Match the position of the end-effector (T_0_3_fun(q)*[0 0 0 1]') in the
% base coordinates to a specific position p3

% p3_val = [0.06 0.04 0.02 1]';
% p3_resp = vpasolve( A_0_3_fun(q)*[0 0 0 1]' == p3_val );
% 
% fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
% vpa(A_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]'))*[0 0 0 1]';
% 
% p3_val = [0.0847 0.0123 -0.005 1]';
% p3_resp = vpasolve( A_0_3_fun(q)*[0 0 0 1]' == p3_val );
% 
% fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
% vpa(A_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]'))*[0 0 0 1]';


