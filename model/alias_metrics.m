function [aliasFull_dB, alias20k_dB, H1] = alias_metrics(y, fs, freq)

    y = y(:)';                                  
    N = length(y);

    y_hann = hann(N)' .* y;
    fftPoints = 10*N;


    Y = fft(y_hann, fftPoints) / fftPoints;
    df = fs / fftPoints;

    idx1 = round(freq/df) + 1;                              % Posicion del fundamental en Y
    H1 = 2 * abs(Y(idx1)) * fftPoints / sum(hann(N));

    numArm = floor((fs/2) / freq);                          % cuántos armónicos caben por debajo de 24 kHz (22)

    for i= 1:numArm

        idx = round(i*freq/df) + 1;
        Y(idx-100 : idx+100) = 0;                           % Borramos los valores de los bins cercanos a harmonicos para aplanar el pico
    
    end
    
    idx20k = round(20000/df) + 1;                        % bin de 20 kHz (misma fórmula de siempre)

    aliasFull = sum(abs(Y(1:fftPoints/2)).^2);            % de 0 a fs/2
    alias20k  = sum(abs(Y(1:idx20k)).^2);                 % de 0 a 20 kHz

    aliasFull_dB = 10*log10(aliasFull);
    alias20k_dB  = 10*log10(alias20k);

end