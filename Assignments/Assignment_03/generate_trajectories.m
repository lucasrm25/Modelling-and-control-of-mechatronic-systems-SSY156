

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

[tq,qv,dqv,ddqv] = traj_planning(qs, ddqc,dts);


figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(tq,qv, 'LineWidth',2); 
grid on, xlabel 't', ylabel 'q'
legend({'q1','q2','q3'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(tq(1:end-1),dqv , 'LineWidth',2 );
grid on, xlabel 't', ylabel 'dq'
ax3 = subplot(3,1,3);
plot(tq(1:end-2), ddqv , 'LineWidth',2);
grid on, xlabel 't', ylabel 'ddq'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Configuration space trajectory planning')
% fp.savefig('trajectory_planning')



%% Question 6 - Operational Space Trajectory Planning


pes = [0   0.098 -0.078 -0.1;
       2   0.098 -0.078 -0.1;
       7   0.2 0 0;
       27   0.2 0 0];

ddpec = 0.5;
dts = 0.1;
[te,pev,dpev,ddpev] = traj_planning(pes, ddpec, dts);

figure('Color','white','Position',[228   414   839   394]);
ax1 = subplot(3,1,1);
plot(te,pev, 'LineWidth',2); 
grid on, xlabel 't', ylabel 'x_e'
legend({'x','y','z'},'Location','southeast')
ax2 = subplot(3,1,2);
plot(te(1:end-1),dpev , 'LineWidth',2 );
grid on, xlabel 't', ylabel 'dx_e'
ax3 = subplot(3,1,3);
plot(te(1:end-2), ddpev , 'LineWidth',2);
grid on, xlabel 't', ylabel 'ddx_e'
linkaxes([ax1,ax2,ax3],'x')
sgtitle('Operational space trajectory planning')
% fp.savefig('trajectory_planning_op_space')


%% Save mat file

save('trajectories_mat.mat', 'qs', 'tq','qv','dqv','ddqv',  'pes', 'te','pev','dpev','ddpev');