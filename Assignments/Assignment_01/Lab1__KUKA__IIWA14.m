clc;clear;

%% Define parameters and variables:
d1=0.36;
d3=0.42;
d5=0.4;
d7=0.131;

q0 = [pi/6; -pi/3; 0; -pi/3; -pi/6; pi/2; pi/2];

q = sym('q',[7,1]);
% v = 
% dq= 


%% Define the transformation matrices and rotation matrices:

DH_table = [0 -pi/2 d1 q(1);
            0  pi/2 0  q(2);
            0  pi/2 d3 q(3);
            0 -pi/2 0  q(4);
            0 -pi/2 d5 q(5);
            0  pi/2 0  q(6);
            0     0 d7 q(7)];

T = rb.T_DH_table(DH_table);    % T{i} = T^{i-1}_i


%% Define Z vector and P vectors:

% Calculate Geometric Jacobian for revolute joints
J_geom = rb.J_geom_rev(T);

digits(10);
J=dre(vpa(J_geom));

matlabFunction(J,'File','Jacobian_KUKA7','Vars',{q},'Optimize',false);
fprintf("Geometric Jacobian written to file!! READY!\n")


%% Define Trajectory

syms t
p= [-0.2;
    -0.0008*t^3 + 0.012*t^2 - 0.2;
    1];
v = diff(p,t);
v_fun = matlabFunction(v);


%% Question 4)c

p_desired = [-0.2 0.2 1];




