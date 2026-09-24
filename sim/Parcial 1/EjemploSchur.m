clear; clc;

P = sdpvar(2,2);
A = [0 1; -1 -2];
LMI1 = A'*P + P*A;
LMI2 = P;

cons = [LMI1 <= 0 , LMI1 >=0];
solvesdp(cons)
Pe = double(P)
eig(Pe)

eig(A'*Pe + Pe*A')