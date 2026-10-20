clear; clc;
iden = 1;
epsilon = 1;

Po = sdpvar(3,3)
F = sdpvar(3,1)

A = [0 1 0;
    0 0 1;
    0 0 0];

C = [1 0 0];

E = [0;
    0;
    1];

alpha = 10;

s11 = Po*A + A'*Po - F*C - C'*F' + 2*alpha*Po;
s12 = Po*E;
s22 = -(1/epsilon)*iden;

LMI1 = [s11 s12; s12' s22];
LMI2 = Po;

cons = [LMI1 <= 0; LMI2 >= 0];
solvesdp(cons);
po = double(Po)
eig(po)
fo = double(F)
Ko = inv(po)*fo

vpa(Ko, 7)

