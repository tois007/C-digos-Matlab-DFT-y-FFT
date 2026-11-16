%% Parcial - Parte 2.1
clear; clc; close all
pareja = 'Perez Torres - Alvarez';                      % apellidos de la pareja
fprintf('%s | %s\n', pareja, string(datetime('now')))   % imprime pareja y fecha/hora

t0 = 0.83e-3; T = 16.69e-3; c = 0;     % los MISMOS valores de la Parte 1 (s, s, mV)
dt = T/10;
tn = t0 + (0:9)*dt;                    % instantes de los 10 puntos de control
xn = [353, 354, 366, 388, 386, -342, -353, -371, -383, -381] - c;   % mV (tabla 1.2)

% Coeficientes de la Parte 1 (DFT directa, seccion 1.4), tal como quedaron en la plantilla
X_mano = [17, 636-2300.8i, 12.8+9.7i, 741-515.5i, 6.2-7.2i, ...
          721, 6.2+7.2i, 741+515.5i, 12.8-9.7i, 663+2300.8i];
X_dft = midft(xn);     % DFT programada en clase
X_fft = fft(xn);       % FFT de MATLAB
% Tabla: k | Re mano | Im mano | Re fft | Im fft | error absoluto
disp('   k     Re mano    Im mano     Re fft     Im fft    |error|')
disp([(0:9)', real(X_mano)', imag(X_mano)', real(X_fft)', imag(X_fft)', abs(X_mano - X_fft)'])
fprintf('max |midft - fft| = %.3e\n', max(abs(X_dft - X_fft)))   % para 2.1 c

%% Parte 2.2 - FFT 2x5 programada
W10 = exp(-1i*2*pi/10);          % W_10
W5 = W10^2;                      % W_5 = W_10^2
xe = xn(1:2:end);                % pares: x0, x2, x4, x6, x8
xo = xn(2:2:end);                % impares: x1, x3, x5, x7, x9
E = zeros(1,5); O = zeros(1,5);
for k = 0:4
    for r = 0:4
        E(k+1) = E(k+1) + xe(r+1)*W5^(k*r);   % DFT de 5 puntos de los pares
        O(k+1) = O(k+1) + xo(r+1)*W5^(k*r);   % DFT de 5 puntos de los impares
    end
end
X25 = zeros(1,10);
for k = 0:4
    X25(k+1) = E(k+1) + W10^k*O(k+1);          % X[k]   = E[k] + W^k O[k]
    X25(k+6) = E(k+1) - W10^k*O(k+1);          % X[k+5] = E[k] - W^k O[k]
end
E, O                                            % comparar con tablas de 1.5
max(abs(X25 - fft(xn)))                         % coincide con la de MATLAB?

%% Parte 2.3 - Los 10 puntos representan la senal?
N = 10;
tau = linspace(0, T, 500);                  % 500 instantes dentro de la ventana
xr = real(X_fft(1))*ones(size(tau));        % termino k = 0
for k = 1:4                                 % armonicos k = 1..4 (cada uno con su espejo)
    xr = xr + 2*( real(X_fft(k+1))*cos(2*pi*k*tau/T) - imag(X_fft(k+1))*sin(2*pi*k*tau/T) );
end
xr = xr + real(X_fft(6))*cos(2*pi*5*tau/T); % termino k = 5
xr = xr/N + c;                              % se divide entre N y se devuelve c

figure; hold on; grid on
plot(t0 + tau, xr, 'b', 'LineWidth', 1.5)             % reconstruccion (azul)
plot(tn, xn + c, 'ro', 'MarkerFaceColor', 'r')        % puntos de control (rojo)
% puntos medios calculados a mano en 1.7 (A: n=4-5, B: n=1-2)
plot(t0 + [4.5 1.5]*dt, [21.5 335.2], 'ks', 'MarkerSize', 9, 'LineWidth', 1.5)
legend('Reconstruccion (10 puntos)', 'Puntos de control', 'Puntos medios a mano (1.7)')
xlabel('t (s)'); ylabel('x (mV)')
title([pareja ' - Representan los 10 puntos a la senal?'])
xlim([t0 t0+T]); ylim([-500 500])          % mismos limites que la foto de la grafica

fprintf('|X[k]| ventana original:   '); fprintf('%8.2f ', abs(X_fft(1:6))); fprintf('\n')

%% Parte 2.3 c - Ventana desplazada medio paso (t0 + dt/2)
% CAMBIAR: leer en la grafica los 10 valores en los nuevos instantes tn2
t0b = t0 + dt/2;
tn2 = t0b + (0:9)*dt;
fprintf('Nuevos instantes (ms): '); fprintf('%.2f ', tn2*1e3); fprintf('\n')
xn2 = [361, 355, 379, 385, -77, -364, -359, -387, -384, 85];  % mV, leidos en captura Tek 2.5 ms/div
% (t = 16.69 ms cae fuera de pantalla: por periodicidad x(16.69 ms) = x(0) = 85 mV)
if numel(xn2) ~= 10
    disp('2.3c pendiente: llenar xn2 con los 10 valores leidos en la grafica en tn2');
else
xn2 = xn2 - c;

X2 = fft(xn2);
xr2 = real(X2(1))*ones(size(tau));
for k = 1:4
    xr2 = xr2 + 2*( real(X2(k+1))*cos(2*pi*k*tau/T) - imag(X2(k+1))*sin(2*pi*k*tau/T) );
end
xr2 = xr2 + real(X2(6))*cos(2*pi*5*tau/T);
xr2 = xr2/N + c;

figure; hold on; grid on
plot(t0b + tau, xr2, 'b', 'LineWidth', 1.5)
plot(tn2, xn2 + c, 'ro', 'MarkerFaceColor', 'r')
legend('Reconstruccion (ventana t_0 + \Deltat/2)', 'Puntos de control')
xlabel('t (s)'); ylabel('x (mV)')
title([pareja ' - Ventana desplazada medio paso'])
xlim([t0b t0b+T]); ylim([-500 500])

fprintf('|X[k]| ventana desplazada: '); fprintf('%8.2f ', abs(X2(1:6))); fprintf('\n')
end

%% Parte 2.4 - El costo de hacerlo a mano
Ns = [10 100 500 1000 2000];
for j = 1:numel(Ns)
    x = randn(1, Ns(j));
    tic; midft(x); t1 = toc;                    % DFT directa
    tic; fft(x);   t2 = toc;                    % FFT
    fprintf('N = %5d   midft = %.5f s      fft = %.6f s\n', Ns(j), t1, t2)
end
