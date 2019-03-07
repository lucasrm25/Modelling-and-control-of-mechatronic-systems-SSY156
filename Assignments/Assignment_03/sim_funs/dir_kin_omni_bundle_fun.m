function out1 = dir_kin_omni_bundle_fun(in1)
%DIR_KIN_OMNI_BUNDLE_FUN
%    OUT1 = DIR_KIN_OMNI_BUNDLE_FUN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    07-Mar-2019 00:02:23

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
t2 = cos(q1);
t3 = cos(q2);
t7 = pi./2.0;
t4 = q3-t7;
t5 = sin(q1);
t6 = sin(q2);
t8 = sin(t4);
t9 = cos(t4);
out1 = [t2.*t3.*(1.3e1./1.0e2)+t2.*t3.*t9.*(1.3e1./1.0e2)-t2.*t6.*t8.*(1.3e1./1.0e2);t3.*t5.*(1.3e1./1.0e2)+t3.*t5.*t9.*(1.3e1./1.0e2)-t5.*t6.*t8.*(1.3e1./1.0e2);t6.*(-1.3e1./1.0e2)-t3.*t8.*(1.3e1./1.0e2)-t6.*t9.*(1.3e1./1.0e2)];
