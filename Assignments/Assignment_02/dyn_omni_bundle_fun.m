function out1 = dyn_omni_bundle_fun(g,in2,in3,in4)
%DYN_OMNI_BUNDLE_FUN
%    OUT1 = DYN_OMNI_BUNDLE_FUN(G,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    19-Feb-2019 11:19:24

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
t16 = t15.*8.838774297938315e89;
t17 = t16+1.767573853083426e91;
t18 = 1.0./t17;
t19 = sign(dq2);
t20 = t19.*3.232508085336712e93;
t21 = dq3.^2;
t22 = q3.*3.0;
t23 = t2+t22;
t24 = cos(t23);
t25 = t10.*t24.*8.960603379009699e91;
t26 = cos(q2);
t27 = g.*t26.*3.354018684496807e94;
t28 = q2+t5;
t29 = cos(t28);
t30 = t13.*t21.*1.174942325548171e93;
t31 = dq2.*dq3.*t13.*2.349884651096342e93;
out1 = [((dq1.*8.945718578607364e47+xi1.*4.676805239458889e49+sign(dq1).*8.888031727671975e47-dq1.*dq2.*t4.*3.805449751309579e47-dq1.*dq3.*t4.*1.902724875654789e47+dq1.*dq2.*t7.*1.205594405358826e47+dq1.*dq3.*t7.*1.205594405358826e47-dq1.*dq2.*t11.*1.684842909617204e48-dq1.*dq3.*t13.*1.902724875654789e47).*2.0)./(t9.*3.805449751309579e47-cos(t2).*1.684842909617204e48+cos(t6).*1.205594405358826e47+sin(t3).*3.805449751309579e47+1.762865565072591e48);t18.*(dq2.*-3.903373707040618e94+dq3.*6.784958835943363e93-t8.*4.332544856330053e93+t20+t25+t27+t30+t31-xi2.*2.887951113948987e95+xi3.*2.887951113948987e95+dq3.*t9.*3.26665061078507e93-g.*t29.*1.32234585863836e93-t4.*t10.*6.770771965641824e92+t7.*t10.*1.41420388767013e92-t8.*t9.*2.085924269165665e93-t10.*t11.*5.343415884116128e93+t10.*t12.*1.41420388767013e92+t10.*t13.*5.874711627740854e92+t12.*t14.*2.828407775340261e92+t13.*t14.*1.174942325548171e93+t9.*xi3.*1.390417760579857e95).*(-1.0./3.2e2);(t18.*(dq2.*-3.903373707040618e94+dq3.*2.329716495749878e94-t8.*1.487643693119772e94+t20+t25+t27+t30+t31-xi2.*2.887951113948987e95+xi3.*9.916209533127894e95-dq2.*t9.*1.879297783897868e94+dq3.*t9.*6.533301221570139e93-g.*t29.*1.32234585863836e93+t4.*t10.*2.004883925276162e93-t7.*t10.*9.058779012314099e92-t8.*t9.*4.17184853833133e93-t10.*t11.*5.201995495349115e93+t10.*t12.*2.828407775340261e92+t10.*t13.*2.017169583169411e93+t12.*t14.*5.656815550680522e92+t13.*t14.*4.034339166338822e93+t9.*t19.*1.556306348594771e93+t12.*t21.*2.828407775340261e92-t9.*xi2.*1.390417760579857e95+t9.*xi3.*2.780835521159713e95+g.*sin(q2+q3).*2.112405578518306e94-g.*sin(q2-q3).*7.75571617167647e93-t10.*cos(-q3+t2).*1.252262701445018e93+dq2.*dq3.*t12.*5.656815550680522e92))./3.2e2];
