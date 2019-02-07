clc;clear;

%% Define parameters and variables:
d1=0.36;
d3=0.42;
d5=0.4;
d7=0.131;

q = sym('q',[7,1]);

%% Define the transformation matrices and rotation matrices:

DH_table = [0 -pi/2 d1 q(1);
            0  pi/2 0  q(2);
            0  pi/2 d3 q(3);
            0 -pi/2 0  q(4);
            0 -pi/2 d5 q(5);
            0  pi/2 0  q(6);
            0     0 d7 q(7)];

% DH_table (Nx4)  Link_i= [a_i,alpha_i,d_i,theta_i]
T = rb.T_DH_table(DH_table);    % T{i} = T^{i-1}_i


%% Define Analytical Jacobian

T_0_e = eye(4);
for i=1:numel(T)
    T_0_e = T_0_e * T{i};
end

% Calculate Geometric Jacobian for revolute joints
J_geom = rb.J_geom_rev(T);

digits(10);
J=dre(vpa(J_geom));
matlabFunction(J,'File','Jacobian_KUKA7','Vars',{q},'Optimize',false);
matlabFunction(T_0_e,'File','KUKA_FKin','Vars',{q},'Optimize',false);
fprintf("Geometric Jacobian written to file!! READY!\n")


%% Set initial conditions

q0 = [pi/6; -pi/3; 0; -pi/3; -pi/6; pi/2; pi/2];

% Question 4c
% q0 = [pi; -pi/3; 0; -pi/3; -pi/6; pi/2; pi/2];

%% Define controller gain

K = eye(6)*0.1;


