function out1 = dyn_omni_bundle_fun(g,in2,in3,in4)
%DYN_OMNI_BUNDLE_FUN
%    OUT1 = DYN_OMNI_BUNDLE_FUN(G,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    21-Feb-2019 21:33:40

dq1 = in3(1,:);
dq2 = in3(2,:);
dq3 = in3(3,:);
q2 = in2(2,:);
q3 = in2(3,:);
xi1 = in4(1,:);
xi2 = in4(2,:);
xi3 = in4(3,:);
t2 = q2.*2.0;
t3 = q3+t2;
t4 = cos(t3);
t5 = q3.*2.0;
t6 = t2+t5;
t7 = sin(t6);
t8 = sign(dq3);
t9 = sin(q3);
t10 = dq1.^2;
t11 = sin(t2);
t12 = sin(t5);
t13 = cos(q3);
t14 = dq2.^2;
t15 = cos(t5);
t16 = sign(dq2);
t17 = dq3.^2;
t18 = q3.*3.0;
t19 = t2+t18;
t20 = cos(t19);
t21 = cos(q2);
t22 = q2+t5;
t23 = cos(t22);
out1 = [((dq1.*3.670737216836427e47-xi1.*7.307508186654515e47-sign(dq1).*1.249272622281644e44+dq1.*dq2.*t4.*2.541906225031903e46+dq1.*dq3.*t4.*1.270953112515952e46+dq1.*dq2.*t7.*3.660659956588252e46+dq1.*dq3.*t7.*3.660659956588252e46-dq1.*dq2.*t11.*1.151283164794982e46+dq1.*dq3.*t13.*1.270953112515952e46).*(-2.0./3.0))./(t9.*8.473020750106344e45+cos(t2).*3.837610549316606e45-cos(t6).*1.220219985529417e46+sin(t3).*8.473020750106344e45+7.778343658881095e46);-(dq2.*6.004939686887598e92-dq3.*5.991601140880533e92-t8.*9.385068133716634e88-t16.*1.224956487367394e89-xi2.*1.200655850930182e93+xi3.*1.200655850930182e93-dq3.*t9.*9.655659793127635e91-g.*t21.*1.382961454899274e92-g.*t23.*1.305951025892415e91-t4.*t10.*8.017969091078025e90+t7.*t10.*8.413134449031172e89-t8.*t9.*1.512434204209622e88+t10.*t11.*8.616731846191844e90+t10.*t12.*8.413134449031172e89+t10.*t13.*1.044116032320737e91+t12.*t14.*1.682626889806234e90+t13.*t14.*2.088232064641474e91-t10.*t20.*2.423191232129347e90+t13.*t17.*2.088232064641474e91+t9.*xi3.*1.93489589053424e92+dq2.*dq3.*t13.*4.176464129282949e91)./(t15.*1.682626889806234e90+2.196855941234287e92);(dq2.*4.803951749510078e93-dq3.*1.298185829279809e94-t8.*2.033440172592318e90-t16.*9.799651898939149e89-xi2.*9.605246807441457e93+xi3.*2.601432196286666e94+dq2.*t9.*7.741724234509505e92-dq3.*t9.*1.544905566900422e93-g.*t21.*1.106369163919419e93-g.*t23.*1.044760820713932e92+t4.*t10.*7.245629893463965e91+t7.*t10.*4.110022377793067e92-t8.*t9.*2.419894726735395e89+t10.*t11.*7.566436232875969e91+t10.*t12.*1.346101511844987e91+t10.*t13.*2.262261130892703e92-t9.*t16.*1.579245724179294e89+t12.*t14.*2.692203023689975e91+t13.*t14.*4.524522261785406e92+t12.*t17.*1.346101511844987e91-t10.*t20.*1.938552985703478e91+t13.*t17.*1.67058565171318e92-t9.*xi2.*1.547916712427392e93+t9.*xi3.*3.095833424854784e93+g.*sin(q2+q3).*2.117483436422934e93+g.*sin(q2-q3).*9.756582157682673e91+t10.*cos(-q3+t2).*6.096778840347489e90+dq2.*dq3.*t12.*2.692203023689975e91+dq2.*dq3.*t13.*3.341171303426359e92)./(t15.*1.346101511844987e91+1.757484752987429e93)];
