%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This Source file plots the comparison between one of the theoretical modulations and their respective experimental   %
%  modulation                                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utilized for the Final Exam of Communication Systems (EET-50) in the Aeronautics Institute of Technology, ITA, Brazil %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ### ITEM A ### --- Obtaining the Error Probability %%
format long;

% Message %
clear all;
clc;
close all;

% Iteration Parameters % 
EbN0min = -2;
EbN0max = 16;
quantity = 100;
errorCounter = 0;
probError2Pam = [];
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantity:EbN0max;

% 2-PAM Modulation:
m = 2;
Es = 1;                 % Symbol Energy Es = (1^2 + 1^2)/2 = 1;
N = 1000;               % Monte Carlo Iterations
messageLength = 1280;   % Length of the transmitted message (symbols)
message = randi([0, m-1], [1, messageLength]);  % Message

% Let's plot the behavior of the modulation as Eb/N0 increases.

for q = 1:quantity
    % Total Error Counter obtained after Monte Carlo %
    errorCounter = 0; 
    %%% Monte Carlo Iterative Loop %%%
    for i = 1:N
        
        % Select Eb/N0 for constructed noise %
        EbN0 = EbN0vector(q);                               
        pamMessage = pammod(message, m);     % PAM Modulation  

        % Construct Symbol Energy-based Noise and Add it %
        noise = (sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(message));
        noisyPAMMessage = pamMessage + noise;

        % PAM Demodulation %
        demodulatedPAMMessage = pamdemod(noisyPAMMessage, m);             
    
        % Analysis: Differences between the received signal and the original signal.
        errorCounter = errorCounter + sum(demodulatedPAMMessage ~= message);
    end
    % Error Probability Computation
    probError2Pam(q) = errorCounter/(messageLength*N);
end

% Theoretical 2-PAM Plot %%

% Basic Parameters
m = 2; % m-PAM
% Data Set for Theoretical Plot % 

EbN0linear = 10.^(EbN0vector./10);

% Theoretical Probability 
argument = sqrt(6*(log2(m)/(m^2-1)).*(EbN0linear));
Pe_t = 2*(1 - 1/m)*qfunc(argument); % Theoretical Probability Formula

% Figure Plot
figure();
semilogy(EbN0vector, Pe_t, 'b');
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Error Probability');
title('2-PAM Error Probability');
semilogy(EbN0vector(1:length(probError2Pam)), probError2Pam, 'r');
axis([-2, 10, 10^-5, 1]);
legend('Theoretical', 'Experimental');
hold off;



%% ### ITEM B ### --- 4-PAM Modulation Calculation %%

% Iteration Parameters % 
EbN0min = -2;
EbN0max = 16;
quantity = 100;
errorCounter = 0;
probError2Pam = [];
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantity:EbN0max;
probError4Pam = [];
% 4-PAM Modulation:
m = 4;
message = randi([0, m-1], [1, messageLength]);
% Symbol Energy: 3^2 + 1^2 + 1^2 + 3^2 = 20 --> Es = 20/4.
Es = 5;  % Symbol Energy

% Let's plot the behavior of the modulation as Eb/N0 increases.
for q = 1:quantity
    errorCounter = 0;
    for i = 1:N
                                                                       
        EbN0 = EbN0vector(q);                                                              
        
        % Message Normalization is performed in the noise. %
        pamMessage = pammod(message, m);
        noise = (sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(message));
        noisyPAMMessage = pamMessage + noise;
        demodulatedPAMMessage = pamdemod(noisyPAMMessage, m);
    
        % Analysis: Differences between the received signal and the original signal.
        errorCounter = errorCounter + sum(demodulatedPAMMessage ~= message);     
    end
    probError4Pam(q) = errorCounter/(messageLength*N);
end

% Theoretical 4-PAM Plot %

% Basic Parameters
m = 4; % m-PAM
EbN0linear = 10.^(EbN0vector./10);

% Calculate the theoretical error probability
argument = sqrt(6*(log2(m)/(m^2-1)).*(EbN0linear));

% Symbol Energy Normalization.
%Eb = Es/K => eS = KEb --> M-pam: 
Pe_t = 2*(1 - 1/m)*qfunc(argument); % Q Function

% Theoretical curve plotted
figure();
semilogy(EbN0vector, Pe_t, 'black');
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Error Probability');
title('4-PAM Error Probability');
semilogy(EbN0vector(1:length(probError4Pam)), probError4Pam, 'green');
axis([-2, 14, 10^-5, 1]);
legend('Theoretical', 'Experimental');
hold off;

%% ### ITEM C ### --- QAM Modulation Error Probability Calculation %%


% Iteration Parameters % 
EbN0min = -2;
EbN0max = 16;
errorCounter = 0;
probError4QAM = [];
% 4-QAM Modulation:
m = 4;
message = randi([0, m-1], [1, messageLength]);
% Calculation of Symbol Energy: It is twice the sqrt(4) - PAM, thus twice
% the 2-PAM
Es = 2;

% Let's plot the behavior of the modulation as Eb/N0 increases.
for q = 1:quantity
    errorCounter = 0;  
    for i = 1:N
        EbN0 = EbN0vector(q); 
 % Modulation, noise, and demodulation %
        qamMessage = qammod(message, m); 
        noise = (sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(message)) ... 
            + 1j*(sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(message));
        noisyQAMMessage = qamMessage + noise; 
        demodulatedQAMMessage = qamdemod(noisyQAMMessage, m);
    
        % Analysis: Differences between the received signal and the original signal.
        errorCounter = errorCounter + sum(demodulatedQAMMessage ~= message);

    end
    probError4QAM(q) = errorCounter/(messageLength*N);
end


% 4-QAM Curve Plot %

EbN0linear = 10.^(EbN0vector./10);

% Calculate the theoretical error probability
argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_t = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2; % Q Function

% Theoretical curve plotted
figure();
semilogy(EbN0vector, Pe_t);
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Error Probability');
title('4-QAM Error Probability');
semilogy(EbN0vector(1:length(probError4QAM)), probError4QAM, 'blue');
axis([-2, 10, 10^-5, 1]);
legend('Theoretical', 'Experimental');
hold off;

%% ### ITEM D ### - FSK Modulation %%
m = 2;         % 2-FSK Modulation
k = log2(m);   % Bits
EbNo = 5;      % Eb/No (dB)
Fs = 16;       % Sampling Rate
nsamp = 8;     % Number of Samples per Symbol
freqsep = 10;  % Frequency Separation

messageLength = 128;
N = 100;
message = randi([0, m-1], [1, messageLength]);

% Iteration Parameters % 
EbN0min = -2;
EbN0max = 16;
probError2FSK = [];
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantity:EbN0max;
Es = 1; % The symbol energy of the BFSK signal is the symbol amplitude itself.

% 2FSK Modulation
for q = 1:quantity
    errorRate = 0;
    for i = 1:N
    
        EbN0 = EbN0vector(q);                                                                 
        mFSK = fskmod(message, m, freqsep, nsamp, Fs);
        mFSKNoisy = awgn(mFSK, (EbN0) + 10*log10(k) - 10*log10(nsamp), 'measured', [], 'dB');
        mFSKDemod = fskdemod(mFSKNoisy, m, freqsep, nsamp, Fs);

        % For 2FSK modulation, the existence of 1 bit makes BER = SER %
        [num, BER] = biterr(message, mFSKDemod);

        errorRate = errorRate + BER;
        
    end
    % Experimental Error Probability %
    probError2FSK(q) = errorRate/N;

    % Theoretical Error Probability Vector (in 2-FSK, BER = SER)%
    BER_theoretical = berawgn(EbN0, 'fsk', m, 'noncoherent');
    Pe_t(q) = BER_theoretical;
end

% Plotting the Curve for 2-FSK
figure();
semilogy(EbN0vector(1:length(Pe_t)), Pe_t);
hold on;
semilogy(EbN0vector(1:length(probError2FSK)), probError2FSK);
grid on;
xlabel('Eb/No (dB)');
ylabel('Error Probability');
title('2-FSK Error Probability');
legend('Theoretical', 'Experimental');
hold off;