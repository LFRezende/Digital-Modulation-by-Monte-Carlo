%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This source file modulates a message into different modulation techniques, aiming to maximing the transmition rate of %
%   data, whilst maintaining a threshold for Bit Error Rate for the message, utilizing Monte Carlo's Algorithm         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;
clear all;
clc;
close all;

% Construction of the Initial and Common Data
quantity = 100;                % Points on the Graph
EbN0min = -2;                    % Minimum Eb/N0
EbN0max = 22;                    % Maximum Eb/N0
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantity:EbN0max; % Construction of the Eb/N0 Vector
EbN0linear = 10.^(EbN0vector./10);% Obtaining its linear equivalent


%% Theoretical Curve for BPSK Modulation %%
% The BPSK constellation is similar (scaled) to the 2PAM constellation.
% Therefore, its error probability has the same formula as 2-PAM.

m = 2;
argument = sqrt(6*(log2(m)/(m^2-1)).*(EbN0linear));
Pe_tBPSK = 2*(1 - 1/m)*qfunc(argument);

%%%% Theoretical Approximation for Comparison - Bit Error Probability %%%%
Pe_tBPSK_bit = Pe_tBPSK/(log2(m));

%% Theoretical Curve for QPSK Modulation %%

% For QPSK modulation, its constellation is identical to the 4QAM constellation.
% Therefore, its error probability has the same formula as 4QAM.

m = 4;
argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_tQPSK = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2;
%%%% Theoretical Approximation for Comparison - Bit Error Probability %%%%
Pe_tQPSK_bit = Pe_tQPSK/(log2(m));

%% Theoretical Curve for 8-PSK Modulation %%

% From M = 8 onwards, the double integral of the error probability calculation
% is no longer separable. But from M = 8 it is already possible to compute using
% the following expression: 2Q(sqrt(2*log2(m)*Eb/N0)*sin(pi/m));

m = 8;
argument = sqrt(2*log2(m).*(EbN0linear))*sin(pi/m);
Pe_t8PSK = 2*qfunc(argument);
%%%% Theoretical Approximation for Comparison - Bit Error Probability %%%%
Pe_t8PSK_bit = Pe_t8PSK/(log2(m));

%% Theoretical Curve for 16-QAM Modulation %%
m = 16;                         % 16-QAM Modulation

% Calculate the theoretical error probability of 16-QAM with complete error
% probability function.

argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_t16 = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2;
%%%% Theoretical Approximation for Comparison - Bit Error Probability %%%%
Pe_t16_bit = Pe_t16/(log2(m));

%% Theoretical Curve for 64-QAM Modulation %%

m = 64;                     % 64-QAM Modulation

% Calculate the theoretical error probability of 64-QAM with complete error
% probability function.

argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_t64 = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2;
%%%% Theoretical Approximation for Comparison - Bit Error Probability %%%%
Pe_t64_bit = Pe_t64/(log2(m));

%% Development of Experimental Modulations %%
% The modulation selection system must be such that:
% If the error probability of the next modulation, for that iteration,
% is below Plim, then we should switch to that modulation.

% Construction of the Message %
len_message = 1280;
message = randi([0, 1], [1, len_message]);
length_message = length(message);

% Iteration Parameters %
N = 500; % Monte Carlo Iteration
counter = 0; % Counting Parameter

% Let's set up some working booleans %
isStart = 1;
switchQPSK = 0;
switch8PSK = 0;
switch16PAM = 0;
switch64PAM = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limiting Probability - >>>> QUESTION 4 <<<<%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the limiting probability of the system below %

probLim_BIT = 10^-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q = 1: quantity
   bitError = 0;
   
    % Modulation Switcher %
    % Here, boolean variables are used to indicate which modulation
    % should be selected.
    % Check for Start
    if (q > 1)
        isStart = 0;
    end
    %   If it has started %
    if (isStart == 0)
        % Check if it switched to the designated modulation, and if it didn't
        % check if the probability of the next modulation is lower than the limiting probability.

        if (switchQPSK == 0)
            if (Pe_tQPSK_bit(q) < probLim_BIT )
                switchQPSK = 1;
                EbN0vector(q)
            end
        end
        if (switch8PSK == 0)
            if (Pe_t8PSK_bit(q)< probLim_BIT)
                switch8PSK = 1;
                EbN0vector(q)
            end
        end
        if (switch16PAM == 0)
            if (Pe_t16_bit(q) < probLim_BIT)
                switch16PAM = 1;
                EbN0vector(q)
            end
        end
        if (switch64PAM == 0)
            if (Pe_t64_bit(q) < probLim_BIT)
                switch64PAM = 1;
                EbN0vector(q)
            end
        end
    end
    

messageLength = len_message;
% Check if we should modulate in BPSK (2-PAM) %
if (mux == [0, 0, 0, 0])
    bitError = 0;
    m = 2;
    Es = 1;
    message = randi([0, m-1], [1, messageLength]);

    % BPSK Modulation (2-PAM) %
    for counter = 1: numSimulations
        % BPSK Modulation (2-PAM) %
        modBPSK = pammod(message, m,[], 'gray'); 
        noise = (sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        noisyBPSK = modBPSK + noise;
        demodBPSK = pamdemod(noisyBPSK, m,[], 'gray'); 
        demodBPSKBinary = dec2bin(demodBPSK);
        messageBinary = dec2bin(message);

        bitError = bitError + sum(sum(demodBPSKBinary ~= messageBinary))/(log2(m)*length(messageBinary));
    % --------------------------------------------------- %
    end
end

% Check if we should modulate in QPSK (4-PAM) %
if (mux == [1, 0, 0, 0])
    bitError = 0;
    m = 4;
    Es = 2;
    message = randi([0, m-1], [1, messageLength]);
    for counter = 1:numSimulations
        % QPSK Modulation (4-QAM) %
        modQPSK = qammod(message, m, 'gray');
        noise = (sqrt(Es/(2*(10.^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) + ...
        1j*(sqrt(Es/(2*(10.^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        noisyQPSK = modQPSK + noise;
        demodQPSK = qamdemod(noisyQPSK, m, 'gray'); 
        demodQPSKBinary = dec2bin(demodQPSK);
        messageBinary = dec2bin(message);

        bitError = bitError + sum(sum(demodQPSKBinary ~= messageBinary))/(log2(m)*length(messageBinary));

    % --------------------------------------------------- %
    end
end

% Check if we should modulate in 8PSK %
if (mux == [1, 1, 0, 0])
    bitError = 0;
    m = 8;
    Es = 1;
    message = randi([0, m-1], [1, messageLength]);
    for counter = 1:numSimulations
        % 8PSK Modulation %
        mod8PSK = pskmod(message, m,[], 'gray'); 
        noise = (sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) ... 
            + 1j*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        noisy8PSK = mod8PSK + noise;
        demod8PSK = pskdemod(noisy8PSK, m,[], 'gray'); 
        demod8PSKBinary = dec2bin(demod8PSK);
        messageBinary = dec2bin(message);

        bitError = bitError + sum(sum(demod8PSKBinary ~= messageBinary))/(log2(m)*length(messageBinary));
    % --------------------------------------------------- %
    end
end

% Check if we should modulate in 16-PAM %
if (mux == [1, 1, 1, 0])
    bitError = 0;
    m = 16;
    Es = 5*2;
    message = randi([0, m-1], [1, messageLength]);
    for counter = 1:numSimulations
        % 16-QAM Modulation %
        mod16QAM = qammod(message, m, 'gray'); 
        noise = (1)*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) ...
            + 1j*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        noisy16QAM = mod16QAM + noise;
        demod16QAM = qamdemod(noisy16QAM, m, 'gray'); 
        demod16QAMBinary = dec2bin(demod16QAM);
        messageBinary = dec2bin(message);

        bitError = bitError + (sum(sum(demod16QAMBinary ~= messageBinary)))/(log2(m)*length(messageBinary));
    % --------------------------------------------------- %
    end
end

% Check if we should modulate in 64-PAM %
if (mux == [1, 1, 1, 1])
    bitError = 0;
    m = 64;
    Es = 42;
    message = randi([0, m-1], [1, messageLength]);
    for counter = 1:numSimulations
        % 64-QAM Modulation %
        mod64QAM = qammod(message, m, 'gray');
        noise = (1)*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) ...
            + 1j*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        noisy64QAM = mod64QAM + noise;
        demod64QAM = qamdemod(noisy64QAM, m, 'gray'); 
        demod64QAMBinary = dec2bin(demod64QAM);
        messageBinary = dec2bin(message);

        bitError = bitError + sum(sum(demod64QAMBinary ~= messageBinary))/(log2(m)*length(messageBinary));
    % --------------------------------------------------- %
    end
end

% Add the count of errors found for this SNR %
BER(q) = bitError/(numSimulations);
end

%% Theoretical Curves Plotted %%
figure();
semilogy(EbN0vector(1: length(BER)), BER, 'red-o');       % Experimental Bit Error Rate
hold on;

% Theoretical Modulation Probabilities %
semilogy(EbN0vector, Pe_tBPSK_bit);
semilogy(EbN0vector, Pe_tQPSK_bit);
semilogy(EbN0vector, Pe_t8PSK_bit);
semilogy(EbN0vector, Pe_t16_bit);
semilogy(EbN0vector, Pe_t64_bit);