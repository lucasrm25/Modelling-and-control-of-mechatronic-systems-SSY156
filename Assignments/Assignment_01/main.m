clear all; close all; clc;

L1 = 0.132;
L2 = 0.132;
syms q1 q2 q3

T_0_1 = rb.Tm_DH(0,-sym(pi)/2,0,q1);
T_1_2 = rb.Tm_DH(L1,0,0,q2);
T_2_3 = rb.Tm_DH(L2,0,0,q3-sym(pi)/2);
T_0_3 = T_0_1 * T_1_2 * T_2_3;


%% Question 1

T_0_3_fun = matlabFunction(T_0_3,'Vars',{[q1 q2 q3]});

q = [0 0 0];
T_0_3_fun(q)

q = [0.67 -0.15 2.7];
T_0_3_fun(q)

q = [-0.73 0.25 1.5];
T_0_3_fun(q)


%% Question 2

p3 = [0.06;0.04;0.02;1];
p3_resp = vpasolve( T_0_3_fun([q1 q2 q3])*[0 0 0 1]' == p3 );
fprintf("q1: %.3f\nq2: %.3f\nq1: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(T_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]))*[0 0 0 1]'

p3 = [0.0847;0.0123;-0.005;1];
p3_resp = vpasolve( T_0_3_fun([q1 q2 q3])*[0 0 0 1]' == p3 );
fprintf("q1: %.3f\nq2: %.3f\nq1: %.3f\n\n", p3_resp.q1, p3_resp.q2, p3_resp.q3 )
vpa(T_0_3_fun([p3_resp.q1, p3_resp.q2, p3_resp.q3]))*[0 0 0 1]'


%% Question 3






