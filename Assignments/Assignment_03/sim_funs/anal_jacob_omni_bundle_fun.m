function out1 = anal_jacob_omni_bundle_fun(in1,in2)
%ANAL_JACOB_OMNI_BUNDLE_FUN
%    OUT1 = ANAL_JACOB_OMNI_BUNDLE_FUN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Mar-2019 08:24:45

phi = in2(1,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
theta = in2(2,:);
t2 = q2+q3;
t3 = cos(t2);
t4 = cos(q1);
t5 = sin(t2);
t6 = t5.*(1.3e1./1.0e2);
t7 = cos(q2);
t8 = t7.*(1.3e1./1.0e2);
t9 = t6+t8;
t10 = sin(q1);
t11 = t3.*(1.3e1./1.0e2);
t12 = sin(q2);
t13 = t11-t12.*(1.3e1./1.0e2);
t14 = phi-q1;
t15 = sin(t14);
t16 = cos(theta);
t17 = sin(theta);
t18 = 1.0./t17;
t19 = cos(t14);
t20 = t15.*t18;
out1 = reshape([-t9.*t10,t4.*t9,0.0,1.0,0.0,0.0,t4.*t13,t10.*t13,-t6-t8,-t15.*t16.*t18,t19,t20,t3.*t4.*(1.3e1./1.0e2),t3.*t10.*(1.3e1./1.0e2),-t6,-t15.*t16.*t18,t19,t20],[6,3]);