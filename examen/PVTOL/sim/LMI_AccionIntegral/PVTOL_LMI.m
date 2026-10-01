clear; clc;

%% Parametros de decaimiento
% Se asigna una especificacion distinta a cada canal para obtener
% nueve ganancias calculadas de forma independiente.
alpha_x = 0.40;
alpha_z = 0.90;
alpha_tau = 3.00;
%% Sistema de un canal
A = [0 1;
     0 0];

B = [0;
     1];

C = [1 0];

%% Sistema aumentado con accion integral
Aa = [A zeros(2,1);
      -C 0];

Ba = [B;
      0];

%% Variables de decision para los tres canales
gammax = sdpvar(3,3);
gammaz = sdpvar(3,3);
gammatau = sdpvar(3,3);

Fx = sdpvar(1,3);
Fz = sdpvar(1,3);
Ftau = sdpvar(1,3);

%% Matrices auxiliares
cerox = zeros(3,3);
ceroz = zeros(3,3);
cerotau = zeros(3,3);

idenx = eye(3);
idenz = eye(3);
identau = eye(3);

%% LMI de x
LMIx = [Aa*gammax + gammax*Aa' ...
    - (Ba*Fx + Fx'*Ba') + alpha_x*gammax, cerox;
    cerox, -idenx];

%% LMI de z
LMIz = [Aa*gammaz + gammaz*Aa' ...
    - (Ba*Fz + Fz'*Ba') + alpha_z*gammaz, ceroz;
    ceroz, -idenz];

%% LMI de theta
LMItau = [Aa*gammatau + gammatau*Aa' ...
    - (Ba*Ftau + Ftau'*Ba') + alpha_tau*gammatau, cerotau;
    cerotau, -identau];

%% Restricciones
consx = [LMIx <= 0, gammax >= 1e-6*eye(3)];
consz = [LMIz <= 0, gammaz >= 1e-6*eye(3)];
constau = [LMItau <= 0, gammatau >= 1e-6*eye(3)];

%% Solucion
diagnostics_x = solvesdp(consx);
if diagnostics_x.problem ~= 0
    error('No se pudo resolver la LMI del canal x: %s', ...
        yalmiperror(diagnostics_x.problem));
end

diagnostics_z = solvesdp(consz);
if diagnostics_z.problem ~= 0
    error('No se pudo resolver la LMI del canal z: %s', ...
        yalmiperror(diagnostics_z.problem));
end

diagnostics_tau = solvesdp(constau);
if diagnostics_tau.problem ~= 0
    error('No se pudo resolver la LMI del canal theta: %s', ...
        yalmiperror(diagnostics_tau.problem));
end

%% Recuperacion de ganancias
fx = double(Fx);
gx = double(gammax);
Px = inv(gx);
Kx = fx*Px;

fz = double(Fz);
gz = double(gammaz);
Pz = inv(gz);
Kz = fz*Pz;

ftau = double(Ftau);
gtau = double(gammatau);
Ptau = inv(gtau);
Ktau = ftau*Ptau;

%% Mostrar resultados
disp('Kx =');
disp(Kx);

disp('Kz =');
disp(Kz);

disp('Ktau =');
disp(Ktau);

%% Verificacion de estabilidad
disp('Polos canal x =');
disp(eig(Aa-Ba*Kx));

disp('Polos canal z =');
disp(eig(Aa-Ba*Kz));

disp('Polos canal theta =');
disp(eig(Aa-Ba*Ktau));