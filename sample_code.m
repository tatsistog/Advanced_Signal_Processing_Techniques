%% Information about team
% Name: Omada fwtia
% Avvakoumidou Anna                 AEM: 8888
% Panagiotopoulos Apostolos         AEM: 8888
% Togkousidis Anastasios            AEM: 8920
%% Initialization
clc; clear;

%% Part 1
sig = rdsamp('data\100.dat');
[~, tms] = ann2rr('data\100', 'atr');

plot(sig(1:2000));   % plot first 2000 samples of signal
hold on;

RR = tms(tms<=2000); % Find R-R in the first 2000 samples
scatter(RR, sig(RR));