function [tv,qv,dqv,ddqv] = traj_planning(qs, ddqc, dts)
    syms t
    qtcell = {};
    syms qt(t)
    for i=1:size(qs,1)-1
        for j=2:size(qs,2)
            qi = qs(i,j);
            qf = qs(i+1,j);
            dt = qs(i+1,1)-qs(i,1);
            if ddqc < 4*abs(qf-qi)/dt^2
                ddqc_n = 4*abs(qf-qi)/dt^2;
                fprintf('Max acceleration too low, changing to: %.2f\n',ddqc_n);
            else
                ddqc_n = ddqc;
            end
            ddqc_n = ddqc_n * sign(qf-qi);
            tc = dt/2 - 0.5*sqrt( (dt^2 *ddqc_n - 4*(qf-qi))/ddqc_n );
            if isnan(tc), tc=0; end

            qt(t) = piecewise( t<0,                  qi,...
                               (0<=t) & (t<=tc),     qi + 0.5*ddqc_n*t^2, ...
                               (tc<t) & (t<=dt-tc),  qi + ddqc_n*tc*(t-tc/2),...
                               (dt-tc<t) & (t<=dt),  qf - 0.5*ddqc_n*(dt-t)^2,...
                               t>dt,                 qf);
            qtcell{i,j-1} = @(t) double(qt(t-qs(i,1)));
        end
    end

    qv = [];
    tv=0:dts:qs(end,1);
    for i=1:numel(tv)
        idx =  min(find(tv(i) < qs(:,1)))-1;
        if isempty(idx)
            idx = size(qs,1)-1; 
        end
        qv(:,i) = cellfun(@(c) c(tv(i)),qtcell(idx,:),'UniformOutput',true)';
    end

    ddt = @(x,dt) conv2(x,[1 -1]/dt,'valid');

    dqv  =  ddt(qv,dts);
    ddqv = ddt(ddt(qv,dts),dts);
end








