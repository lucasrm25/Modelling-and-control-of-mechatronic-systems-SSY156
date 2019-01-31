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
n_links = size(DH_table,1);

T = rb.T_DH_table(DH_table);    % T{i} = T^{i-1}_i


%% Define Z vector and P vectors:

T_0_i{1} = T{1};
z{1} = [0 0 1]';                        % z0
p{1} = [0 0 0 1]';
for i=2:n_links
    T_0_i{i} = T_0_i{i-1} * T{i};
    z{i} = T_0_i{i-1}(1:3,1:3) * z{1};    % z{2} = z1 = R_0_1 * z0
    p{i} = T_0_i{i-1}          * p{1};    % p{2} = p1 = A_0_1 * p0
end
pe = T{n_links} * p{n_links-1}; 


sel = [eye(3) zeros(3,1)];              % remove last row from vector
J_geom = sym(zeros(6,n_links-1));
for i=1:n_links
    J_geom(1:3,i) = cross(z{i},sel*(pe-p{i}));
    J_geom(4:6,i) = z{i};
end
J_geom = simplify(J_geom);

digits(10);
J=dre(vpa(J_geom));

matlabFunction(J,'File','Jacobian_KUKA7','Vars',{q},'Optimize',false);
fprintf("Geometric Jacobian written to file!! READY!\n")


% matlabFunction(P7,'File','KUKA_FKin','Vars',{q});