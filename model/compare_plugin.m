amp = "0.0981";
ref = load("plugin_ref_amp" + amp + ".mat");
out = load("plugin_out_amp" + amp + ".mat");

uki_K = ref.uki_k;
uki_comp = ref.uki_comp;

yK = double(out.yK);
yA = double(out.yA);

eK = yK(1:end-1) - uki_K(2:end);                        
eA = yA(1:end-2) - uki_comp(3:end);

errK = max(abs(eK));
errK_dB = 20*log10(errK / max(abs(ref.uki_k)));

errA = max(abs(eA));
errA_dB = 20*log10(errA / max(abs(ref.uki_comp)));

fprintf("amp = %s | K: max|e| = %.3g V (%.1f dB) | ADAA1: max|e| = %.3g V (%.1f dB)\n", amp, errK, errK_dB, errA, errA_dB);

%% Figura: plugin frente a MATLAB
% Arriba: las dos salidas superpuestas en los primeros 2 ms (ya alineadas).
% Abajo: el error a lo largo de toda la senal (si creciera con el tiempo,
% indicaria un problema de memoria/estado en el plugin).
fs = ref.fs;
nZoom = round(2e-3 * fs);                 % 2 ms = 96 muestras a 48 kHz

% Senales alineadas (mismo desfase que en eK y eA)
mK = uki_K(2:end);     pK = yK(1:end-1);
mA = uki_comp(3:end);  pA = yA(1:end-2);

tZoom = (0:nZoom-1) / fs * 1e3;           % tiempo en ms
tK = (0:length(eK)-1) / fs;               % tiempo en s
tA = (0:length(eA)-1) / fs;

figure("Name", "Plugin vs MATLAB, amp = " + amp);
tiledlayout(2, 2, "TileSpacing", "compact");

nexttile;
plot(tZoom, mK(1:nZoom), "-", "Color", [0 0.45 0.74]); hold on;
plot(tZoom, pK(1:nZoom), "o", "Color", [0.85 0.33 0.10], "MarkerSize", 4);
hold off; grid on;
title("K-method"); xlabel("t (ms)"); ylabel("V");
legend("MATLAB", "Plugin", "Location", "best");

nexttile;
plot(tZoom, mA(1:nZoom), "-", "Color", [0 0.45 0.74]); hold on;
plot(tZoom, pA(1:nZoom), "o", "Color", [0.85 0.33 0.10], "MarkerSize", 4);
hold off; grid on;
title("ADAA1"); xlabel("t (ms)"); ylabel("V");
legend("MATLAB", "Plugin", "Location", "best");

nexttile;
plot(tK, eK, "Color", [0 0.45 0.74]); grid on;
title(sprintf("Error K: max %.3g V (%.1f dB)", errK, errK_dB));
xlabel("t (s)"); ylabel("Plugin - MATLAB (V)");

nexttile;
plot(tA, eA, "Color", [0 0.45 0.74]); grid on;
title(sprintf("Error ADAA1: max %.3g V (%.1f dB)", errA, errA_dB));
xlabel("t (s)"); ylabel("Plugin - MATLAB (V)");


%% Métricas del estándar: MATLAB frente a plugin
freq = 1073;

[aFull_mK, a20k_mK, H1_mK] = alias_metrics(uki_K, fs, freq);
[aFull_pK, a20k_pK, H1_pK] = alias_metrics(yK, fs, freq);
[aFull_mA, a20k_mA, H1_mA] = alias_metrics(uki_comp, fs, freq);
[aFull_pA, a20k_pA, H1_pA] = alias_metrics(yA, fs, freq);

fprintf("\namp = %s\t| alias full \t| alias 20k \t| H1\n", amp);

fprintf("K MATLAB  \t| %8.2f \t| %8.2f \t| %.4f V\n", aFull_mK, a20k_mK, H1_mK);
fprintf("K Plugin  \t| %8.2f \t| %8.2f \t| %.4f V\n", aFull_pK, a20k_pK, H1_pK);
fprintf("ADAA1 MATLAB  \t| %8.2f \t| %8.2f \t| %.4f V\n", aFull_mA, a20k_mA, H1_mA);
fprintf("ADAA1 Plugin  \t| %8.2f \t| %8.2f \t| %.4f V\n", aFull_pA, a20k_pA, H1_pA);

