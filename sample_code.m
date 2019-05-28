%% Information about team
% Name: Omada fwtia
% Avvakoumidou Anna                 AEM: 8888
% Panagiotopoulos Apostolos         AEM: 8888
% Togkousidis Anastasios            AEM: 8920

%% Initialization
clc; clear;

% Array of data
data = [];
data(1).pathdat = 'data\113.dat';
data(1).pathhea = 'data\113.hea';
data(1).pathatr = 'data\113';

data(2).pathdat = 'data\116.dat';
data(2).pathhea = 'data\116.hea';
data(2).pathatr = 'data\116';

data(3).pathdat = 'data\203.dat';
data(3).pathhea = 'data\203.hea';
data(3).pathatr = 'data\203';

data(4).pathdat = 'data\231.dat';
data(4).pathhea = 'data\231.hea';
data(4).pathatr = 'data\231';

data(5).pathdat = 'data\234.dat';
data(5).pathhea = 'data\234.hea';
data(5).pathatr = 'data\234';

%% Part 1 - Algorithms evalutation
signal_struct = [];
frequencies = [];
RR_struct = [];
spect_struct = [];

for i = 1:5
    % ---------- Signals over time  --------------------------------------%    
    
    [sig, fs, tm] = rdsamp(data(i).pathdat);
    [rr, tms] = ann2rr(data(i).pathatr, 'atr');
    
    % Saving up values
    signal_struct(i).ecg1 = sig(:,1);
    signal_struct(i).ecg2 = sig(:,2);
    frequencies(i) = fs;
    RR_struct(i).rr = rr;
    
    % Plot ECG
    figure;
    plot(sig(1:2000,:)); hold on;       % plot first 2000 samples of signal
    scatter(tms(tms<=2000), sig(tms(tms<=2000)));
    mytitle = strcat("ECG from data ", string(i)) ;
    xlabel('Samples');
    title(mytitle);
    
    % ---------- Short Time Fourier Transform ----------------------------%
    % Parameters
    Nx = size(sig,1);
    nsc = floor(Nx/5000);
    nov = floor(nsc*0.75);
    nff = max(256,2^nextpow2(nsc));
    
    % Plot spectrogramm
    figure;
    s1 = spectrogram(sig(:,1),gausswin(nsc),nov,nff,fs);
    spectrogram(sig(:,1),gausswin(nsc),nov,nff,fs);
    mytitle = strcat("Spectrogramm from data ", string(i)) ;
    title(mytitle);
    
    % Saving up spectrogramm matrix
    spect_struct(i).s = s1;
    
    % -------------- Wigner Distribution Function ------------------------%
    samps = 10000;
    W = zeros(samps);
    for j = 0: (650000/samps - 1)
        W = W + mywigner(sig((j*samps+1):(j+1)*samps,1));
    end
    W = W / j;
    
    % Plotting  WDF
    figure;
    h = surf(W,'edgecolor', 'none');
    view(2); colorbar;
    mytitle = strcat("Wigner Distribution Function from data ", string(i)) ;
    title(mytitle);
    
    % End of first loop
    fprintf('Plots of data %d. \n',i);
    fprintf('Press any key to continue...\n');
    pause;
    close all;
    
end
%% Part 2 - STFT - Parameters must get fixed

% Short Time Fourier Transform
Nx = size(sig,1);
nsc = floor(Nx/5000);
nov = floor(nsc/8);
nff = max(256,2^nextpow2(nsc));

figure;
s1 = spectrogram(sig(:,1),gausswin(nsc),nov,nff,fs);
spectrogram(sig(:,1),gausswin(nsc),nov,nff,fs);

%% Part 3 - WDF - Parameters must get fixed

% Wigner Distribution Function
W = mywigner(sig(1:1000,1));

% Plotting
figure;
h = surf(W);
set(h,'LineStyle','none');
view(2); colorbar;

%% Part 4 - Wavelet transform
wt = cwt(sig(1:1000,1), 'morse', fs);
cwt(sig(1:400,1), 'bump', fs);

%% Part 5 - Hilbert-Huang Transform
