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

p3 = [0.06 0.04 0.02 1]';
p3_resp = vpasolve( A_0_3_fun(q)*[0 0 0 1]' == p3 );

fprintf("q1: %.3f\nq2: %.3f\nq3: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(A_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]'))*[0 0 0 1]';

p3 = [0.0847 0.0123 -0.005 1]';
p3_resp = vpasolve( A_0_3_fun(q)*[0 0 0 1]' == p3 );

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

pl{1} = p1;
pl{2} = p2 + A{1} * A{2} * [-L1/2 0 0 1].';
pl{3} = p3 + A{1} * A{2} * A{3} * [-L2/2 0 0 1].';

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

% Calculate generalized forces
xi = sym('xi',[3,1]);

L = T - U;

fp.m2latex(pl{1})
fp.m2latex(pl{2})
fp.m2latex(pl{3})
fp.m2latex(J_l{1})
fp.m2latex(J_l{2})
fp.m2latex(J_l{3})
fp.m2latex(A_0_{1})
fp.m2latex(A_0_{2})
fp.m2latex(A_0_{3})
fp.m2latex(T)
fp.m2latex(U)


      
      
      
      

