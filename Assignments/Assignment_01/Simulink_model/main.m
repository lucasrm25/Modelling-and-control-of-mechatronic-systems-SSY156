clear all; close all; clc;
%% OmniBundle parameters

L1 = 0.132;
L2 = 0.132;
syms q1 q2 q3
q = [q1 q2 q3];

T_0_1 = rb.T_DH(0,-sym(pi)/2,0,q1);
T_1_2 = rb.T_DH(L1,0,0,q2);
T_2_3 = rb.T_DH(L2,0,0,q3-sym(pi)/2);
T_0_3 = T_0_1 * T_1_2 * T_2_3;


%% Question 1

T_0_3_fun = matlabFunction(T_0_3,'Vars',{[q1 q2 q3]});

q_val = [0 0 0];
T_0_3_fun(q_val)*[0 0 0 1]'

q_val = [0.67 -0.15 2.7];
T_0_3_fun(q_val)*[0 0 0 1]'

q_val = [-0.73 0.25 1.5];
T_0_3_fun(q_val)*[0 0 0 1]'


%% Question 2

% Match the position of the end-effector (T_0_3_fun(q)*[0 0 0 1]') in the
% base coordinates to a specific position p3

p3 = [0.06 0.04 0.02 1]';
p3_resp = vpasolve( T_0_3_fun(q)*[0 0 0 1]' == p3 );

fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(T_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]))*[0 0 0 1]'

p3 = [0.0847 0.0123 -0.005 1]';
p3_resp = vpasolve( T_0_3_fun(q)*[0 0 0 1]' == p3 );

fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(T_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]))*[0 0 0 1]'


%% Question 3

% 3 revolute joints
z0 = [0 0 1]';
z1 = T_0_1(1:3,1:3)*z0;
z2 = T_0_1(1:3,1:3)*T_1_2(1:3,1:3)*z0;

p0 = [0 0 0 1]';            % selection of 4th column
p1 = T_0_1 * p0;
p2 = T_0_1 * T_1_2 * p0;
pe = T_0_1 * T_1_2 * T_2_3 * p0;

sel = [eye(3) zeros(3,1)];  % remove last row from vector

J_geom = [cross(z0,sel*(pe-p0)), cross(z1,sel*(pe-p1)), cross(z2,sel*(pe-p2));
          z0, z1, z2];
J_geom = simplify(J_geom)

syms phi theta psi
eulZYZ = [phi theta psi].';     % rotm2eul(T_0_3(1:3,1:3),'ZYZ')

J_analit = rb.T_zyz2geom(eulZYZ) \ J_geom;
J_analit = simplify(J_analit)

