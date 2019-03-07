
%% Animate

% config
delay_factor = 2;       % 1 = no delay, 0.5 = half of real speed
fps = 60;

ta = 0:1/fps:simdata.time(end);
data_anim = [ta' interp1(simdata.time, simdata.signals(1).values, ta') ];


% figure('Color','white','Position',[491   446   611   239])
% plot(data_anim(:,1), data_anim(:,2:4)), grid on;
% xlabel 'time [s]', ylabel 'Joint angles q [rad]';
% legend({'q_1','q_2','q_3'})
% fp.savefig('sim')


% Animation

% data_anim = data.q;
    
% close all;
p1_fun = [0 0 0]';
p2_fun = matlabFunction( subs(p2(1:3),unkownvars, unkownvars_val), 'Vars', {q});
p3_fun = matlabFunction( subs(pe(1:3),unkownvars, unkownvars_val), 'Vars', {q});

figure('Color','white'), hold on, grid on, axis equal;

endeff_pos = p3_fun(qs(:,2:4)');
scatter3( endeff_pos(1,:), endeff_pos(2,:), endeff_pos(3,:), 'red' ,'filled');

lin2 = plot3( [0 [1 0 0]*p2_fun(q0)] ,[0 [0 1 0]*p2_fun(q0)], [0 [0 0 1]*p2_fun(q0)] ,'LineWidth',3);
lin3 = plot3( [0 [1 0 0]*p3_fun(q0)] ,[0 [0 1 0]*p3_fun(q0)], [0 [0 0 1]*p3_fun(q0)] ,'LineWidth',3);
xlim([0 2*0.132])
ylim([-0.15 0.15])
zlim([-0.15 0.1])
view(45,15)
tic
for i=2:length(data_anim)
    p2_val = p2_fun(data_anim(i,2:4)');
    p3_val = p3_fun(data_anim(i,2:4)');
    lin2.XData = [0 p2_val(1)];
    lin2.YData = [0 p2_val(2)];
    lin2.ZData = [0 p2_val(3)];
    lin3.XData = [p2_val(1) p3_val(1)];
    lin3.YData = [p2_val(2) p3_val(2)];
    lin3.ZData = [p2_val(3) p3_val(3)];
    
    plot3( [[1 0 0]*p3_fun(data_anim(i-1,2:4)') [1 0 0]*p3_fun(data_anim(i,2:4)')], ...
           [[0 1 0]*p3_fun(data_anim(i-1,2:4)') [0 1 0]*p3_fun(data_anim(i,2:4)')], ...
           [[0 0 1]*p3_fun(data_anim(i-1,2:4)') [0 0 1]*p3_fun(data_anim(i,2:4)')],...
           'LineWidth',3, 'Color', 'black');
    
    drawnow();
    pause(data_anim(i,1)-toc*delay_factor)
    if mod(i,10)==0
        fprintf('%f\n',data_anim(i,1))
    end
end
