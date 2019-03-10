clear all; close all; clc;



syms Tm Tv km Kp Kv  
syms s

G = Kp * Kv * (1 + Tv*s)/s * km/(1 + Tm*s) *1/s
H = 1 + 1/Kp*s

T = simplify(G / (1 + G*H) );
pretty(T)

Tpzc = simplify(subs(T,Tv,Tm));
pretty(Tpzc)




T = simplify(1 / (1 + G*H) );
pretty(T)

Tpzc = simplify(subs(T,Tv,Tm));
pretty(Tpzc)







km = 2;
Tm = 7.2;

Kp = linspace(0,10,50);
Kv = linspace(0.01,0.3,50);
[KP,KV] = meshgrid(Kp,Kv);

DRF = KP.*KV;
TR  = max(Tm, 2./(KV*km));
WN = (KP.*KV.*km).^0.5;
ZETA = (KV .* km).^0.5 ./ (2 * (KP).^0.5);

TAU = 1./(WN.*ZETA);
figure
surf(KP,KV,TAU)
xlabel 'Kp', ylabel 'Kv', zlabel 'TAU'

figure
subplot(2,2,1)
surf(KP,KV,DRF)
xlabel 'Kp', ylabel 'Kv', zlabel 'DRF'
subplot(2,2,2)
surf(KP,KV,TR)
xlabel 'Kp', ylabel 'Kv', zlabel 'ORT'
subplot(2,2,3)
surf(KP,KV,WN)
xlabel 'Kp', ylabel 'Kv', zlabel 'WN'
subplot(2,2,4)
surf(KP,KV,ZETA)
xlabel 'Kp', ylabel 'Kv', zlabel 'ZETA'



% Tm < 2/(KV*km)
% 
% KV < 2/(Tm*km)


% Kp Kv = wn2 / km
% Kv = 2 ksi wn / km

% wn2 = Kp Kv km
% ksi = wn / 2 Kp = sqrt(Kp Kv Km) / 2 Kp = sqrt(Kv km) / (2 Kp sqrt(Kp)


% ksi wn = 2 / Kv.km






