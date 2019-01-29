close all; clear all; clc;

%% User input

theta = 30;
L1 = 1;
L2 = 1;
L3 = 1;
L4 = 0.2;
L5 = 1;
L6 = 1;
L7 = 0.2;


%% Matrices
Tm = @(r,p) [r p; zeros(1,3) 1];



syms L3 L2 L4 L7 L6 L5 theta
rz = @(th) [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
rx = @(th) [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
ry = @(th) [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];

A_O_AUX = Tm(rz(theta),[-L3;L2;L6]) * Tm(ry(pi/2),[L4;0;L7]) * Tm(rz(pi/2),[0;0;0])
vpa(A_O_AUX)

A_O_AUX = Tm(rz(theta),[-L3;L2;L6]) * Tm(rx(pi/2),[0;-L5/2;0])
vpa(A_O_AUX)






A_O_AUX =  Tm(rotz(theta),[-L3;L2;L6]);
A_O_LA  =  A_O_AUX * Tm(rotx(-90),[0;+L5/2;0]);
A_O_RA  =  A_O_AUX * Tm(rotx(+90),[0;-L5/2;0]);
A_O_C =  A_O_AUX * Tm(roty(pi/2),[L4;0;L7]) * Tm(rotz(pi/2),[0;0;0]);


%% Plot

% plt = @(V) quiver3( repmat(V(1,1),3,1), repmat(V(2,1),3,1), repmat(V(3,1),3,1),...
%                     V(1,2:4)'-repmat(V(1,1),3,1), V(2,2:4)'-repmat(V(2,1),3,1), V(3,2:4)'-repmat(V(3,1),3,1),...
%                     'LineWidth',3, 'MaxHeadSize',0.5);
                
figure('Color','white'); grid on, hold on, axis square;
xlabel 'x', ylabel 'y', zlabel 'z';
axs = [zeros(3,1) 0.2*eye(3) ;ones(1,4)];
plt(axs)
plt(A_O_AUX*axs)
plt(A_O_LA*axs)
plt(A_O_RA*axs)
plt(A_O_C*axs)


%% Help Functions

function plt(V)
    quiver3( repmat(V(1,1),3,1), repmat(V(2,1),3,1), repmat(V(3,1),3,1),...
                    V(1,2:4)'-repmat(V(1,1),3,1), V(2,2:4)'-repmat(V(2,1),3,1), V(3,2:4)'-repmat(V(3,1),3,1),...
                    'LineWidth',3, 'MaxHeadSize',0.5);
    text(V(1,2),V(2,2),V(3,2),'x');
    text(V(1,3),V(2,3),V(3,3),'y');
    text(V(1,4),V(2,4),V(3,4),'z');
end