function out1 = inv_dyn_omni_bundle_fun(g,in2,in3,in4)
%INV_DYN_OMNI_BUNDLE_FUN
%    OUT1 = INV_DYN_OMNI_BUNDLE_FUN(G,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Mar-2019 08:24:43

ddq1 = in4(1,:);
ddq2 = in4(2,:);
ddq3 = in4(3,:);
dq1 = in3(1,:);
dq2 = in3(2,:);
dq3 = in3(3,:);
q2 = in2(2,:);
q3 = in2(3,:);
t2 = q2.*2.0;
t3 = q3+t2;
t4 = cos(t3);
t5 = q3.*2.0;
t6 = t2+t5;
t7 = sin(t6);
t8 = dq1.^2;
t9 = sin(t2);
t10 = sin(q3);
t11 = cos(q3);
t12 = ddq3.*1.238e-3;
t13 = q2+q3;
t14 = sin(t13);
out1 = [ddq1.*7.477375e-3+dq1.*8.9e-3+ddq1.*t10.*6.76e-4+ddq1.*cos(t2).*9.58375e-4-ddq1.*cos(t6).*2.19e-4+ddq1.*sin(t3).*6.76e-4+dq1.*dq2.*t4.*1.352e-3+dq1.*dq3.*t4.*6.76e-4+dq1.*dq2.*t7.*4.38e-4+dq1.*dq3.*t7.*4.38e-4-dq1.*dq2.*t9.*1.91675e-3+dq1.*dq3.*t11.*6.76e-4;ddq2.*5.25475e-3+dq2.*(1.7e1./1.0e3)+t12+ddq2.*t10.*1.352e-3+ddq3.*t10.*6.76e-4-g.*t14.*5.2e-3-t4.*t8.*6.76e-4-t7.*t8.*2.19e-4+t8.*t9.*9.58375e-4+dq3.^2.*t11.*6.76e-4-g.*cos(q2).*1.755e-2+dq2.*dq3.*t11.*1.352e-3;ddq2.*1.238e-3+dq3.*5.8e-3+t12+ddq2.*t10.*6.76e-4-g.*t14.*5.2e-3-t4.*t8.*3.38e-4-t7.*t8.*2.19e-4-t8.*t11.*3.38e-4-dq2.^2.*t11.*6.76e-4];
