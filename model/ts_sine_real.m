%--------------------------------------------------------------------------
% Project: Tube Screamer Clipping Stage Simulation
% Method:  K-method with First-Order Antiderivative Antialiasing (ADAA)
%
% Description:
%   This script simulates the clipping stage of an Ibanez Tube Screamer
%   pedal using the K-method for state-space discretization. It compares
%   the performance of the standard K-method against the proposed 
%   ADAA-enhanced version to demonstrate aliasing reduction capabilities.
%
%   Comparison modes:
%   1. Sine wave input: Analysis of harmonic content and noise floor.
%   2. Chirp input: Spectrogram visualization of aliasing artifacts.
%
% Requirements:
%   - binarySearchInterp.m (Must be in the same folder or path)
%   - char.txt (Optional: if spiceModel = 1)
%
%--------------------------------------------------------------------------

clear; close all;

% binarySearchInterp.m esta en el repo de Jedla (submodulo jedla_ref)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'jedla_ref'));

%% Configuration
% Set to 1 for chirp input (spectrogram visualization)
% Set to 0 for sine input (spectrum analysis & harmonic removal)
isChirp = 0; 

%% Plot Defaults
% Store old defaults to restore later if needed
oldDefaultFontSize   = get(groot, 'DefaultAxesFontSize');
oldDefaultGridLine   = get(groot, 'DefaultAxesXGrid');
oldDefaultLineColor  = get(groot, 'DefaultLineColor');
oldDefaultLineWidth  = get(groot, 'DefaultLineLineWidth');

% Set local defaults for publication-quality plots
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultAxesXGrid', 'on', 'DefaultAxesYGrid', 'on');
set(groot, 'DefaultAxesColorOrder', [0 0 0]);   % RGB triplet for black
set(groot, 'DefaultLineLineWidth', 1.5);

%% Generating / Loading Parameters
spiceModel = 0; % 1: Load from file, 0: Generate ideal diode char

% Upsampling parameters
upsample_ideal = 8;  % High oversampling for reference (8x or 16x)
upsample_k = 1;      % Comparison target for standard method
upsample_new = 1;    % ADAA method (usually lower to prove efficiency)

fs = 48e3;
R1 = 4700;        % 4.7k (R6)

paramSet = "real";
if(paramSet=="real")
    Rgain = 200e3;    % 200K (R7)
    C1 = 47e-9;       % 47nF
    N = 1.752;        %diode ideality factor (1N4148, LTspice model)
    UT = 25.85e-3;    % 25.85mV

elseif(paramSet=="jedla") % ideal diode (Jedla)
    Rgain = 0;        % 0 (R7)
    C1 = 47e-6;       % 47uF
    N = 1;            
    UT = 26e-3;       % 26mV

else
    error("Unknown paramSet: %s (use ""real"" or ""jedla"")", paramSet);
    
end
R2 = 51e3+Rgain;  % 51k (R5) + drive potentiometer (R7)
C2 = 51e-12;      % 51pF
Is0 = 2.52e-9;    % 2.52nA



freq = 1073;      % Test frequency for sine
amp = 0.0981;     % Input amplitude (change between 0.9798 and 0.0981)
dur = 2;          % Simulation duration in seconds

% State-space matrices
% w' = Aw + Bu + Cy
% x  = Dw + Eu + Fy
% y  = f(x)
A = [-1/(R1*C1) 0; -1/(R1*C2) -1/(R2*C2)];
B = [1/(R1*C1); 1/(R1*C2)];
C = [0; -1/C2];
D = [0 1];

% Initial parameters for calculations (Standard)
h = 2 * fs * upsample_new;
I2 = eye(2);
invA_init = inv(h*I2 - A);
G = invA_init * C;
K = D * G;

xMax = 1;
pointsPerVoltage = 10000;

if spiceModel == 1
    spModel = load("char.txt");
    y = spModel(:,2);
    x = spModel(:,1);
else
    x = 0 : 1/pointsPerVoltage : xMax;
    y = Is0 * (exp(x/(N*UT)) - exp(-x/(N*UT)));
end

% Pre-calculate lookup tables
p = x - K * y;          % p domain for K-method
Y = cumtrapz(p, y);     % Antiderivative for ADAA

% --- Debug Plots for Nonlinear Functions ---
figure(1);
plot(x, y);
xlim([0, max(x)]);
xlabel("x (Voltage)", 'FontSize', 16); ylabel("y (Current)", 'FontSize', 16);
title("Original Nonlinear Function f(x)");

figure(2);
plot(p, y);
xlim([0, max(p)]);
xlabel("p (Pseudo-variable)", 'FontSize', 16); ylabel("y (Current)", 'FontSize', 16);
title("Transformed Nonlinear Function (K-method)");

figure(3);
plot(p, Y);
xlabel("p (Pseudo-variable)", 'FontSize', 16); ylabel("Y (Antiderivative)", 'FontSize', 16);
xlim([0, max(p)]);
title("Antiderivative F(p) for ADAA");


%% 1. Modified K-method with ADAA
fprintf('Running Modified K-method (ADAA)...\n');

% Compensation for ADAA frequency response
% Reescalado ADAA1 (tesis de Jedla, ecs. 4.5-4.9): h~ = 4/(3T) -> fs~ = 2/3*fs.
% Antes era fs_ = fs (sin reescalar, como en jedla_ref/ts.m); ver registro 2026-10-21 §3.
fs_ = fs*2/3;
h_ = 2 * fs_ * upsample_new;

% Pre-calculated matrices
invA = inv(h_ * I2 - A);
G_tilde = invA * C;
H_tilde = invA * (h_ * I2 + A);
J_tilde = invA * B;
K_tilde = D * G_tilde;

% Recalculate p vector for the specific step size
p_vec = x - K_tilde * y;
Y_vec = cumtrapz(p_vec, y);

% --- Parameters and Pre-allocation ---
N_new = floor(dur * fs * upsample_new);
t_new = (0 : N_new-1) / (fs * upsample_new);

% Generate Input
if isChirp
    u_raw = amp * chirp(t_new, 20, dur, fs/2, "linear");
else
    u_raw = amp * sin(2*pi*freq*t_new);
end

u_padded = [0, 0, u_raw(1:end-1)];
if length(u_padded) < N_new
    u_padded(N_new) = 0; 
end

% Pre-allocation
w = zeros(2, N_new);
pk = zeros(2, N_new);
pn = zeros(1, N_new);
yn = zeros(1, N_new);
uki_comp = zeros(1, N_new);
tolerance = eps;

% --- Simulation Loop (ADAA) ---
for i = 4 : N_new
    % Calculate pk (predictor)
    pk(:, i) = H_tilde * w(:, i-1)  + ...
               J_tilde * (u_padded(i) + 0.5*(u_padded(i-1) + u_padded(i-2))) + ...
               G_tilde * yn(i-1);
    
    pn(i) = D * pk(:, i);
    pnBuffer = pn(i-1);
    delta_pn = pn(i) - pnBuffer;
    
    % ADAA Logic
    if abs(delta_pn) < tolerance
        % Fallback to standard lookup if change is too small (avoid div by zero)
        midpoint_p = (pn(i) + pnBuffer) * 0.5;
        yi = binarySearchInterp({p_vec}, {y}, abs(midpoint_p));
        if midpoint_p < 0
             yi = -yi; % Odd function symmetry
        end
    else
        % ADAA step: [F(p_n) - F(p_n-1)] / (p_n - p_n-1)
        Fpn       = binarySearchInterp({p_vec}, {Y_vec}, abs(pn(i)));
        FpnBuffer = binarySearchInterp({p_vec}, {Y_vec}, abs(pnBuffer));
        yi = (Fpn - FpnBuffer) / delta_pn;
    end
    
    yn(i) = yi;
    
    % State update
    w(:, i) = G_tilde * yn(i) + 0.5 * (pk(:, i) + pk(:, i-1));
    
    % Output calculation
    uki_comp(i) = w(2, i) + u_padded(i);
end

% Downsampling (if needed)
if upsample_new > 1
    uki_comp = resample(uki_comp, 1, upsample_new);
end


%% 2. Standard K-method (Reference)
fprintf('Running Standard K-method...\n');

fs_k = fs * upsample_k;
h = 2 * fs_k;

% Pre-calculated matrices
invA = inv(h * I2 - A);
G = invA * C;
H = invA * (h * I2 + A);
J = invA * B;
K = D * G;

p_k_vec = x - K * y;

% --- Time and Input ---
N_k = floor(dur * fs_k);
t_k = (0 : N_k-1) / fs_k;

if isChirp
    u_k_raw = amp * chirp(t_k, 20, dur, fs/2, "linear");
else
    u_k_raw = amp * sin(2*pi*freq*t_k);
end

u_k = [0, u_k_raw(1:end-1)];
if length(u_k) < N_k
    u_k(N_k) = 0;
end

% Pre-allocation
pk = zeros(2, N_k);
pn = zeros(1, N_k);
w_k = zeros(2, N_k);
yn_k = zeros(1, N_k);
uki_k = zeros(1, N_k);

% --- Simulation Loop (Standard) ---
for i = 2 : N_k
    pk(:, i) = H * w_k(:, i-1) + ...
               J * (u_k(i) + u_k(i-1)) + ...
               G * yn_k(i-1);
    
    pn(i) = D * pk(:, i);
    
    % Standard Lookup
    yi = binarySearchInterp({p_k_vec}, {y}, abs(pn(i)));
    if pn(i) < 0
        yi = -yi;
    end
    
    yn_k(i) = yi;
    w_k(:, i) = G * yn_k(i) + pk(:, i);
    uki_k(i) = w_k(2, i) + u_k(i);
end

% Downsampling
if upsample_k > 1
    uki_k = resample(uki_k, 1, upsample_k, 50);
end


%% 3. Ideal Case (High Oversampling)
fprintf('Running Ideal Case (Reference)...\n');

fs_ideal = fs * upsample_ideal;
h = 2 * fs_ideal;

% Pre-calculated matrices
invA = inv(h * I2 - A);
G = invA * C;
H = invA * (h * I2 + A);
J = invA * B;
K = D * G;

p_ideal_vec = x - K * y;

% --- Time and Input ---
N_ideal = floor(dur * fs_ideal);
t_ideal = (0 : N_ideal-1) / fs_ideal;

if isChirp
    u_ideal_raw = amp * chirp(t_ideal, 20, dur, fs/2, "linear");
else
    u_ideal_raw = amp * sin(2*pi*freq*t_ideal);
end

u_ideal = [0, u_ideal_raw(1:end-1)];
if length(u_ideal) < N_ideal
    u_ideal(N_ideal) = 0;
end

% Pre-allocation
pk = zeros(2, N_ideal);
pn = zeros(1, N_ideal);
w_ideal = zeros(2, N_ideal);
yn = zeros(1, N_ideal);
uki_ideal = zeros(1, N_ideal);

% --- Simulation Loop (Ideal) ---
for i = 2 : N_ideal
    pk(:, i) = H * w_ideal(:, i-1) + ...
               J * (u_ideal(i) + u_ideal(i-1)) + ...
               G * yn(i-1);
    
    pn(i) = D * pk(:, i);
    
    yi = binarySearchInterp({p_ideal_vec}, {y}, abs(pn(i)));
    if pn(i) < 0
        yi = -yi;
    end
    
    yn(i) = yi;
    w_ideal(:, i) = G * yn(i) + pk(:, i);
    uki_ideal(i) = w_ideal(2, i) + u_ideal(i);
end

% Downsampling
if upsample_ideal > 1
    uki_ideal = resample(uki_ideal, 1, upsample_ideal, 50); %cambio frente a Jedla, pues el usaba n standard (n=10), n=50 reduce error en ideal.
end


%% General Plots (Time Domain)
t = 0 : 1/fs : dur - 1/fs;

% Ensure output vectors are same length as t
len_t = length(t);
uki_comp = uki_comp(1:len_t);
uki_k = uki_k(1:len_t);
uki_ideal = uki_ideal(1:len_t);

figure(5);
plot(t, uki_comp, 'Color', 'cyan', 'LineWidth', 1); hold on;
plot(t, uki_k, 'Color', 'g', 'LineWidth', 1);
plot(t, uki_ideal, 'Color', 'b', 'LineWidth', 1);

legend("Modified K-method (ADAA)", "Standard K-method", "Ideal", 'Location', 'best', 'FontSize', 16);
xlabel('Time [s]', 'FontSize', 16); ylabel('Amplitude [V]', 'FontSize', 16);
title('Output Signals (Time Domain)');


%% Calculation of Comparison Numbers
power_ideal = sum(uki_ideal.^2);
power_new = sum(uki_comp.^2);
power_k = sum(uki_k.^2);

fprintf('\nComparison Results:\n');
fprintf('Ideal Power: %.4f\n', power_ideal);
fprintf('New Power: %.4f\n', power_new);
fprintf('K Power: %.4f\n', power_k);


%% Conditional Plots based on Input Type
if isChirp == 0
    % --- SINE INPUT PLOTS ---
    fprintf('Generating Spectrum Plots (Sine Input)...\n');
    
    fftPoints = length(uki_comp) * 10;
    idx20k = round(20000 / (fs/fftPoints)) + 1;
    f = 0 : 1 / fftPoints * fs : fs - 1 / fftPoints * fs;
    
    % 1. Full Spectrum Plot (Figure 8)
    figure(8);
    plot(f, 20*log10(abs(fft(hann(length(uki_comp))'.*uki_comp, fftPoints)/ fftPoints) + eps), 'Color', 'cyan', 'LineWidth', 1.5); hold on;
    plot(f, 20*log10(abs(fft(hann(length(uki_k))'.*uki_k, fftPoints)/ fftPoints) + eps), 'Color', 'r', 'LineWidth', 1);
    plot(f, 20*log10(abs(fft(hann(length(uki_ideal))'.*uki_ideal, fftPoints)/ fftPoints) + eps), 'Color', 'g', 'LineWidth', 1);
    
    xlim([0 fs/2]);
    xlabel("Frequency [Hz]", 'FontSize', 16);
    ylabel("Amplitude [dB]", 'FontSize', 16);
    legend("Modified K-method (ADAA)", "Standard K-method", "Ideal", 'Location', 'best', 'FontSize', 16);
    grid on;
    title("Full Frequency Spectrum");

    % 2. Harmonic Removal Calculation
    fft_ideal = fft(hann(length(uki_ideal))'.*uki_ideal, fftPoints);
    fft_new = fft(hann(length(uki_comp))'.*uki_comp, fftPoints);
    fft_k = fft(hann(length(uki_k))'.*uki_k, fftPoints);
    
    fft_abs_new_without_harmonics = abs(fft_new) - abs(fft_ideal);
    fft_abs_k_without_harmonics = abs(fft_k) - abs(fft_ideal);
    
    power_new_without_harmonics = sum(fft_abs_new_without_harmonics.^2);
    power_k_without_harmonics = sum(fft_abs_k_without_harmonics.^2);
    
    fprintf('\nPower without Harmonics (Difference from Ideal):\n');
    fprintf('New: %.4f\n', power_new_without_harmonics);
    fprintf('K: %.4f\n', power_k_without_harmonics);
    
    % 3. Advanced Harmonic Removal for Visualization (Noise Floor Plot)
    uki_comp_fft = fft(hann(length(uki_comp))'.*uki_comp, fftPoints) / fftPoints;
    uki_k_fft = fft(hann(length(uki_k))'.*uki_k, fftPoints)/ fftPoints;
    
    uki_ideal_fft=fft(hann(length(uki_ideal))'.*uki_ideal, fftPoints) / fftPoints;


    point_around_peak = 100;
    bin_freq = fs / fftPoints;
    peak_num = floor((fs/2) / freq);
    
    for i = 1 : peak_num
        center_freq = freq * i;
        center_idx = round(center_freq / bin_freq) + 1;
        
        start_idx = max(1, center_idx - point_around_peak);
        end_idx = min(fftPoints, center_idx + point_around_peak);
        
        if end_idx > start_idx
            nullable_items = start_idx : end_idx;
            uki_comp_fft(nullable_items) = 0;
            uki_k_fft(nullable_items) = 0;
            uki_ideal_fft(nullable_items)=0;
        end
    end
    % This loop removes the fundamental and its harmonics from both
        % spectra to isolate the remaining noise floor.
        %
        % For each harmonic i:
        %   1. Compute its frequency: freq * i.
        %   2. Convert that frequency to an FFT-bin index.
        %   3. Select a window of bins around the harmonic peak.
        %   4. Set those bins to zero in both spectra.
        %
        % The window size is controlled by point_around_peak, while the
        % loop stops at the highest harmonic below the Nyquist frequency

    pot_ADAA = sum(abs(uki_comp_fft(1:fftPoints/2)).^2);
    pot_k = sum(abs(uki_k_fft(1:fftPoints/2)).^2);
    pot_ideal = sum(abs(uki_ideal_fft(1:fftPoints/2)).^2);

    pot_ADAAx20k = sum(abs(uki_comp_fft(1:idx20k)).^2);
    pot_kx20k = sum(abs(uki_k_fft(1:idx20k)).^2);
    pot_idealx20k = sum(abs(uki_ideal_fft(1:idx20k)).^2);

    cum_ADAA = cumsum(abs(uki_comp_fft(1:fftPoints/2)).^2); % cumsum returns a vector: in each position,
    cum_k = cumsum(abs(uki_k_fft(1:fftPoints/2)).^2);       % the sum form 0Hz tot that freq

    pot_cum_ADAAdB = 10*log10(cum_ADAA);
    pot_cum_kdB = 10*log10(cum_k);

    pot_ADAAdB = 10*log10(pot_ADAA);
    pot_kdB = 10*log10(pot_k);
    pot_idealdB = 10*log10(pot_ideal);

    pot_ADAAdBx20k = 10*log10(pot_ADAAx20k);
    pot_kdBx20k = 10*log10(pot_kx20k);
    pot_idealdBx20k = 10*log10(pot_idealx20k);
    
    fprintf('\nAliasing power (full band):\n');
    fprintf('ADAA: %.2f\n', pot_ADAAdB);
    fprintf('K: %.2f\n', pot_kdB);
    fprintf("Ideal: %.2f\n", pot_idealdB);

    fprintf('\nAliasing power (0-20 kHz):\n');
    fprintf('ADAA: %.2f\n', pot_ADAAdBx20k);
    fprintf('K: %.2f\n', pot_kdBx20k);
    fprintf("Ideal: %.2f\n", pot_idealdBx20k);

    figure(10);
    plot(f, 20*log10(abs(uki_comp_fft) + eps), 'Color', 'cyan', 'LineWidth', 1); hold on;
    plot(f, 20*log10(abs(uki_k_fft) + eps), 'Color', 'r', 'LineWidth', 1);
    
    xlim([0 fs/2]);
    xlabel("Frequency [Hz]", 'FontSize', 16);
    ylabel("Amplitude [dB]", 'FontSize', 16);
    legend("Modified K-method (Harmonics Removed)", "K-method (Harmonics Removed)", 'Location', 'best', 'FontSize', 16);
    title("Spectrum with Fundamental and Harmonics Removed (Noise Floor)");
    grid on;

    figure(11);
    plot(f(1:fftPoints/2), pot_cum_ADAAdB, 'Color', 'cyan', 'LineWidth', 1); hold on;
    plot(f(1:fftPoints/2), pot_cum_kdB, 'Color', 'r', 'LineWidth', 1);

    xlim([0 fs/2]);
    xlabel("Frequency [Hz]", 'FontSize', 16);
    ylabel("Cumulative power [dB]", 'FontSize', 16);
    legend("ADAA", "K-method", 'Location', 'best', 'FontSize', 16);
    title("Cumulative aliasing power");
    grid on;

else
    % --- CHIRP INPUT PLOTS (Spectrograms) ---
    fprintf('Generating Spectrograms (Chirp Input)...\n');
    
    % Common Settings
    win_size = 1024;       
    overlap = 900;         
    nfft = 4096;           
    dB_range = [-80 0];    
    fig_pos_1 = [100, 500, 1200, 400]; 
    fig_pos_2 = [100, 100, 1200, 400];
    
    % 1. ADAA Spectrogram
    figure('Name', 'ADAA Spectrogram', 'Color', 'w', 'Position', fig_pos_1);
    [s, f_axis, t_axis] = spectrogram(uki_comp, blackmanharris(win_size), overlap, nfft, fs);
    
    S_dB = 20*log10(abs(s));
    S_dB = S_dB - max(S_dB(:)); 
    
    imagesc(t_axis, f_axis/1000, S_dB); 
    axis xy; 
    colormap(flipud(gray)); 
    clim(dB_range);         
    colorbar;
    
    ylabel('Frequency (kHz)', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 16); 
    ylim([0 fs/2000]); 
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 1.5);
    title('Modified K-method with ADAA');
    
    % 2. Standard K-method Spectrogram
    figure('Name', 'Reference Spectrogram', 'Color', 'w', 'Position', fig_pos_2);
    [s_k, f_axis_k, t_axis_k] = spectrogram(uki_k, blackmanharris(win_size), overlap, nfft, fs);
    
    S_k_dB = 20*log10(abs(s_k));
    S_k_dB = S_k_dB - max(S_k_dB(:)); 
    
    imagesc(t_axis_k, f_axis_k/1000, S_k_dB);
    axis xy; 
    colormap(flipud(gray)); 
    clim(dB_range); 
    colorbar;
    
    ylabel('Frequency (kHz)', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 16);
    ylim([0 fs/2000]); 
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 1.5);
    title('Standard K-method');
    
    % 3. Combined Stacked Spectrogram (Shared Axis, Compact)
    figure('Name', 'Combined Spectrogram', 'Color', 'w', 'Position', [100, 100, 1000, 500]);
    
    % Settings: [left, bottom, width, height]
    h_plot = 0.4; 
    gap = 0.03; % Minimal gap
    
    pos_bottom = [0.13, 0.15, 0.80, h_plot]; 
    pos_top    = [0.13, 0.15 + h_plot + gap, 0.80, h_plot];

    % Top Plot (ADAA) - No X Label
    subplot('Position', pos_top);
    imagesc(t_axis, f_axis/1000, S_dB);
    axis xy; colormap(flipud(gray)); clim(dB_range); colorbar;
    ylim([0 fs/2000]);
    
    % Trick for legend
    hold on;
    h_dummy1 = plot(nan, nan, 'LineStyle', 'none'); 
    legend(h_dummy1, 'ADAA', 'Location', 'southeast', 'FontSize', 16, 'Color', 'w', 'EdgeColor', 'k');
    
    % Remove X ticks and individual Y label
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 1.5, 'XTickLabel', []); 
    
    % Bottom Plot (Standard) - With X Label
    subplot('Position', pos_bottom);
    imagesc(t_axis_k, f_axis_k/1000, S_k_dB);
    axis xy; colormap(flipud(gray)); clim(dB_range); colorbar;
    xlabel('Time (s)', 'FontSize', 16); 
    ylim([0 fs/2000]);
    
    % Trick for legend
    hold on;
    h_dummy2 = plot(nan, nan, 'LineStyle', 'none'); 
    legend(h_dummy2, 'Standard', 'Location', 'southeast', 'FontSize', 16, 'Color', 'w', 'EdgeColor', 'k');
    
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 1.5);

    % Shared Y Axis Title
    % Create an invisible axis covering both plots to center the label
    total_height = 2*h_plot + gap;
    h_label = axes('Position', [0.13, 0.15, 0.70, total_height], 'Visible', 'off');
    h_label.YLabel.Visible = 'on';
    ylabel(h_label, 'Frequency (kHz)', 'FontSize', 16, 'Color', 'k');
end

%% Restore old defaults at the end (Best practice)
set(groot, 'DefaultAxesFontSize', oldDefaultFontSize);
set(groot, 'DefaultAxesXGrid', oldDefaultGridLine);
set(groot, 'DefaultLineColor', oldDefaultLineColor);
set(groot, 'DefaultLineLineWidth', oldDefaultLineWidth);

%% 4. Validación frente a LTspice (parte lineal)

freq_vect = [100 1073 20000];
fSweep = logspace(1,5,4001);   % 10 Hz, 10.02 Hz, ... , 100 kHz
HdB = zeros(size(fSweep));

for k = 1:numel(fSweep)
    s = 1j*2*pi*fSweep(k);
    M = s*eye(2)-A;
    X=M\B;
    Hlin=D*X+1;
    HdB(k)=20*log10(abs(Hlin));
    %fprintf("|H(%g)| = %.2f [dB]\n", fSweep(k), HdB(k));
end

[Hmax, idx] = max(HdB);
fprintf("\nHmax = %.2fdB at freq %gHz\n", Hmax, fSweep(idx));

%% 5. Armónicos: Ideal frente a LTspice y ADAA1 frente a Ideal.

df = fs / fftPoints;        % bin size
Hmag = zeros(1, 15);        % a place to store the harmonics
Hmag_new = zeros(1, 15);    % ADAA1 (fft_new)
for n = 1:15
    fharm=freq*n;
    iBin=round(fharm/df)+1;
    Hmag(n)=abs(fft_ideal(iBin));
    Hmag_new(n)=abs(fft_new(iBin));

end

Hnorm_dB=20*log10(Hmag./Hmag(1));
Hnorm_new_dB=20*log10(Hmag_new./Hmag_new(1));
v_step1 = [
    1.000e+00;
    3.537e-03;
    2.274e-01;
    3.030e-03;
    1.189e-01;
    2.670e-03;
    7.495e-02;
    2.364e-03;
    5.153e-02;
    2.089e-03;
    3.718e-02;
    1.841e-03;
    2.763e-02;
    1.614e-03;
    2.093e-02
];

v1_dB = 20*log10(v_step1);

v_step2 = [
    1.000e+00;
    2.424e-03;
    1.176e-01;
    9.635e-04;
    6.626e-02;
    8.483e-04;
    4.516e-02;
    8.178e-04;
    3.379e-02;
    7.885e-04;
    2.672e-02;
    7.630e-04;
    2.190e-02;
    7.400e-04;
    1.841e-02
];

v2_dB = 20*log10(v_step2);

H1_V = 2*Hmag(1)/sum(hann(length(uki_ideal)));
H1_V_new = 2*Hmag_new(1)/sum(hann(length(uki_comp)));

fprintf("\n| H1_ideal = %.4fV | H1_ADAA = %.4fV | H1_dif(dB) = %.4fdB |\n", H1_V, H1_V_new, 20*log10(abs(H1_V_new/H1_V)));

fprintf("\n_________________________________________________________________________________________________________________\n");

if amp==0.9798
    fprintf("\nstep 2\n");
    vdB=v2_dB;

elseif amp==0.0981
    fprintf("\nstep 1\n");
    vdB=v1_dB;

else
    error("No hay datos de LTspice para amp = %g (usa 0.0981 o 0.9798)", amp);

end 

for n = 1:15
    if mod(n, 2) == 1
        fprintf("\n| n = %g \t| f = %g \t| MATLAB dBc = %.2fdB \t| LTSpice dBc = %.2fdB \t| dif = %.2f \t|", n, freq*n, Hnorm_dB(n), vdB(n), Hnorm_dB(n)-vdB(n));
    end
end

fprintf("\n_________________________________________________________________________________________________________________\n");

for n = 1:15
    if mod(n, 2) == 1
        fprintf("\n| n = %g \t| Hnorm(ideal) = %.4fdB \t| Hnorm(ADAA) = %.4fdB \t| Hnorm_diff = %.4fdB \t|", n, Hnorm_dB(n), Hnorm_new_dB(n), Hnorm_new_dB(n)-Hnorm_dB(n));
    end
end

fileName = sprintf("plugin_ref_amp%.4f.mat", amp);
save(fileName, "u_raw", "uki_k", "uki_comp", "fs", "amp", "-v7");
