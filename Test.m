clear all
close all
clc

% Parámetros del sistema
K = 2;      % Ganancia
tau = 5;    % Constante de tiempo

% Modelo en Espacio de Estados
A = -1/tau;
B = K/tau;
C = 1;
D = 0;

% Crear el objeto del sistema
plant = ss(A, B, C, D);

% Parámetros del MPC
Ts = 0.1;  % Tiempo de muestreo
p = 20;    % Horizonte de Predicción
m = 10;    % Horizonte de Control

% Crear el controlador MPC
mpc_controller = mpc(plant, Ts, p, m);

% Tiempo de Simulación
Tfinal = 50;
t = 0:Ts:Tfinal;

% Referencia (setpoint)
r = ones(size(t));
r(t > 20 & t <= 40) = 3;   % Cambia el setpoint a 3 en t = 20s

% Asegurarse de que r es un vector columna
r = r(:);

% Simular el sistema con el controlador MPC
[y, time, u] = sim(mpc_controller, Tfinal, r);

% Asegurarse de que t, y, y r tengan la misma longitud
if length(y) < length(t)
    t = t(1:length(y));
end

if length(r) < length(t)
    r = [r; r(end)*ones(length(t)-length(r), 1)];
end


% Graficar los resultados
figure;

subplot(2, 1, 1);
plot(time, y, 'b', t, r, 'r--');
xlabel('Time (s)');
ylabel('Output');
legend('Output y(t)', 'Reference r(t)');
title('System Output vs Reference');

subplot(2, 1, 2);
plot(time, u, 'g');
xlabel('Time (s)');
ylabel('Control Input');
title('Control Input u(t)');
