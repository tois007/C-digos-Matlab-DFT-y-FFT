# C-digos-Matlab-DFT-y-FFT
En este se encuentran los archivos de cada una de las transformadas en Matlab
function X = midft(x)
% MIDFT  DFT directa por definicion: X[k] = sum_n x[n] * exp(-i*2*pi*k*n/N)
% (usar la version de clase si la tienen; esta es equivalente)
    N = numel(x);
    X = zeros(1, N);
    for k = 0:N-1
        for n = 0:N-1
            X(k+1) = X(k+1) + x(n+1) * exp(-1i*2*pi*k*n/N);
        end
    end
end
