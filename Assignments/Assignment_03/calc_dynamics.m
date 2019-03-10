
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
simplify(jacobian(gradL,q)*dq + jacobian(gradL,dq)*ddq - gradient(L,q) - (B*ddq + h + gvec));

 

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


%% Simulink help functions

dyn_omni_bundle = subs(DYNEQ, unkownvars, unkownvars_val);
[Adyn,bdyn]=equationsToMatrix(dyn_omni_bundle,ddq);
matlabFunction( simplify(Adyn\bdyn), 'Vars', {g, q,dq,xi}, 'File', 'sim_funs/dyn_omni_bundle_fun');

[Adyn,bdyn]=equationsToMatrix(dyn_omni_bundle,xi);
matlabFunction( simplify(Adyn\bdyn), 'Vars', {g, q,dq,ddq}, 'File', 'sim_funs/inv_dyn_omni_bundle_fun');

matlabFunction( subs(pe(1:3),unkownvars, unkownvars_val), 'Vars', {q}, 'File', 'sim_funs/dir_kin_omni_bundle_fun');
matlabFunction( subs(A{1}*A{2}*A{3},unkownvars, unkownvars_val), 'Vars', {q}, 'File', 'sim_funs/A_0_e_omni_bundle_fun');
matlabFunction( subs(J_analit,unkownvars, unkownvars_val), 'Vars', {q,eulZYZ}, 'File', 'sim_funs/anal_jacob_omni_bundle_fun');

dJ_analit_P = sym(zeros(3));
for i=1:9
    dJ_analit_P(i) = jacobian(J_analit(i),q)*dq;
end
matlabFunction( subs( dJ_analit_P, unkownvars, unkownvars_val), 'Vars', {q,dq}, 'File', 'sim_funs/dJ_analit_P_omni_bundle_fun');

addpath('sim_funs/')

