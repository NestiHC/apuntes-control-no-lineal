clear;
clc;

%% =========================================================
%  LESO - Calculo de ganancias mediante LMI
% ==========================================================

%% Matrices del observador

A = [0 1 0;
    0 0 1;
    0 0 0];

C = [1 0 0];

E = [0;
    0;
    1];

%% =========================================================
%  Parametros de diseño
% ==========================================================

alpha = 1.0;
epsilon = 1.0;

%% =========================================================
%  Variables de decision
% ==========================================================

Po = sdpvar(3,3,'symmetric');
Fo = sdpvar(3,1);

%% =========================================================
%  Restricciones LMI
% ==========================================================

LMI1 = Po >= 1e-6*eye(3);

LMI2 = [ ...
    Po*A + A'*Po - Fo*C - C'*Fo' + 2*alpha*Po,   Po*E;
    E'*Po,                                       -epsilon*eye(1) ...
    ] <= -1e-6*eye(4);

%% =========================================================
%  Resolver
% ==========================================================

Constraints = [LMI1, LMI2];

options = sdpsettings('solver','sedumi','verbose',1);

sol = optimize(Constraints,[],options);

%% =========================================================
%  Verificar solucion
% ==========================================================

if sol.problem ~= 0
    disp('No se encontro una solucion factible.');
    disp(sol.info);
    return;
end

%% =========================================================
%  Obtener resultados
% ==========================================================

Po_val = value(Po);
Fo_val = value(Fo);

%% Ganancia del observador

Ko = Po_val\Fo_val;

%% =========================================================
%  Mostrar resultados
% ==========================================================

disp('============================================');
disp('Matriz Po:');
disp(Po_val);

disp('============================================');
disp('Matriz Fo:');
disp(Fo_val);

disp('============================================');
disp('Ganancia del observador Ko:');
disp(Ko);

disp('============================================');

Ko1 = Ko(1);
Ko2 = Ko(2);
Ko3 = Ko(3);

fprintf('Ko1 = %.10f\n',Ko1);
fprintf('Ko2 = %.10f\n',Ko2);
fprintf('Ko3 = %.10f\n',Ko3);