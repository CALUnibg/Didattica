syms A A_u A_e g C rho k q_e W T_e h T real positive
syms f1(h, q_e) f2(h, T, q_e, W, T_e)
f1 = (q_e - A_u * sqrt(2*g*h))/A;
f2 = 1/(C*rho*A*h) * (-C*rho*T*q_e + rho*q_e*(C*T_e+(q_e^2/(2*A_e^2))+g*h) - rho*A_u*sqrt(2*g*h)*g*h + W - k*(T-T_e));
J = jacobian([f1 f2], [h, T, q_e, W, T_e]);
Jeq = simplify(subs(J, [A A_u A_e g C rho k q_e W T_e h T], [0.08 2.5e-5 1e-4 9.8 4184 1e3 10 4.95e-5 1035 20 0.2  24.7672]));
Jeqn = double(Jeq);

s = tf('s');
G11 = Jeqn(1,3) / (s - Jeqn(1,1));
G12 = Jeqn(1,4) / (s - Jeqn(1,1));
G21 = Jeqn(2,3) / (s - Jeqn(2,2));
G22 = Jeqn(2,4) / (s - Jeqn(2,2));
M11 = Jeqn(1,5) / (s - Jeqn(1,1));
M21 = Jeqn(2,5) / (s - Jeqn(2,2));
G = [G11, G12; G21 G22];
M = [M11; M21];

dcgain(G) * [4.95e-5; 1035] + dcgain(M) * 20
