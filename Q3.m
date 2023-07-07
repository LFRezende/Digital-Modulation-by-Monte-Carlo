%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This source file modulates a message into different modulation techniques, aiming to maximing the transmition rate of %
%   data, whilst maintaining a threshold for Symbol Error Rate for the message, utilizing Monte Carlo's Algorithm         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;
clear all;
clc;
close all;
% We will develop theoretical curves for the Symbol Error Rate (SER) of the BPSK, QPSK,
% 8PSK, 16QAM, and 64QAM modulations.

%% Theoretical Probabilities of the Modulations

% Construction of Initial Data Common to All
quantity = 100;                    % Points on the Graph
EbN0min = -2;                      % Minimum Eb/N0
EbN0max = 22;                      % Maximum Eb/N0
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantity:EbN0max; % Construction of the Eb/N0 Vector
EbN0linear = 10.^(EbN0vector./10); % Conversion to linear scale


%% Theoretical Curve of BPSK Modulation %% 

% The BPSK constellation is similar (scaled) to the 2PAM constellation.
% Therefore, its error probability has the same formula as 2-PAM.
% Theoretical probability of 2-PAM == BPSK
m = 2;
argument = sqrt(6*(log2(m)/(m^2-1)).*(EbN0linear));
Pe_tBPSK = 2*(1 - 1/m)*qfunc(argument); 

%% Theoretical Curve of QPSK Modulation %% 

% For QPSK modulation, its constellation is identical to the 4QAM constellation.
% Therefore, its error probability has the same formula as 4QAM.

m = 4;
argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_tQPSK = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2;

%% Theoretical Curve of 8-PSK Modulation %% 

% Starting from M = 8, the double integral for calculating the error probability
% is no longer separable. But from M = 8, it is already possible to compute it using
% the following expression: 2Q(sqrt(2*log2(m)*Eb/N0)*sin(pi/m));
m = 8;
argument = sqrt(2*log2(m).*(EbN0linear))*sin(pi/m);
Pe_t8PSK = 2*qfunc(argument);

%% Theoretical Curve of 16-QAM Modulation %% 
m = 16;                            % 16-QAM Modulation

% Calculate the theoretical error probability of 16-QAM with the complete
% error probability function.

argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_t16 = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2;


%% Theoretical Curve of 64-QAM Modulation %% 

m = 64;                            % 64-QAM Modulation

% Calculate the theoretical error probability of 64-QAM with the complete
% error probability function.

argument = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_t64 = 4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2;


%% Development of Experimental Modulations %%
% The modulation selection decision system should be such that:
% If the error probability of the next modulation, for that iteration,
% is below Plim, then we should switch to that modulation.

% Message Construction %
len_message = 1280;
message = randi([0, 1], [1, len_message]);
length = length(message);

% Iteration Parameters %
N = 1000; % Monte Carlo iteration
counter = 0; % Counter parameter for counting

% Let's set up some working booleans %
isStart = 1;
switchQPSK = 0;
switch8PSK = 0;
switch16QAM = 0;
switch64QAM = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limiting Probability %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the system's limiting probability below %

probLim = 10^-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q = 1: quantity
   errors = 0;
   
    % Modulation Switcher %
    % Here, boolean variables are used to indicate which modulation
    % should be selected.
    % Check if it is the start
    if (q > 1)
        isStart = 0;
    end
    % If it has started %
    if (isStart == 0)
        % Check if it has switched to the designated modulation, and if not
        % check if the probability of the next modulation is
        % lower than the limiting probability. %

        if (switchQPSK == 0)
            if (Pe_tQPSK(q) < probLim )
                switchQPSK = 1;
                EbN0vector(q)
            end
        end
        if (switch8PSK == 0)
            if (Pe_t8PSK(q) < probLim)
                switch8PSK = 1;
                EbN0vector(q)
            end
        end
        if (switch16QAM == 0)
            if (Pe_t16(q) < probLim)
                switch16QAM = 1;
                EbN0vector(q)
            end
        end
        if (switch64QAM == 0)
            if (Pe_t64(q) < probLim)
                switch64QAM = 1;
                EbN0vector(q)
            end
        end
    end
    % Storing the Switching Variables in a Vector %
    mux = [switchQPSK, switch8PSK, switch16QAM, switch64QAM];
    
    % Check if we should modulate with BPSK (2-PAM) %
    if (mux == [0, 0, 0, 0])
        errors = 0;
        m = 2;
        Es = 1;
        message = randi([0, m-1], [1, len_message]);

        % % BPSK Modulation === 2-PAM % %
        for counter = 1: N
            
            % BPSK Modulation == 2-PAM %

            % According to item 1 - modulate, add noise, demodulate, accumulate error %
            mBPSK = pammod(message, m); 
            noise = (sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
            mBPSKnoisy = mBPSK + noise;
            mBPSKdemod = pamdemod(mBPSKnoisy, m); 
            errors = errors + sum(mBPSKdemod ~= message);
        % --------------------------------------------------- %
        end
    end
    % Check if we should modulate with QPSK (4-PAM) %
    if (mux == [1, 0, 0, 0])
        errors = 0;
        m = 4;
        Es = 2;
        message = randi([0, m-1], [1, len_message]);
        for counter = 1:N

    % % QPSK Modulation === 4-QAM % %
        % According to item 1 - modulate, add noise, demodulate, accumulate error %
        m4QAM = qammod(message, m);                % QPSK Modulation == 4-QAM
        noise = (sqrt(Es/(2*(10.^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) + ...
        1j*(sqrt(Es/(2*(10.^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        m4QAMnoisy = m4QAM + noise;
        m4QAMdemod = qamdemod(m4QAMnoisy, m); 
        errors = errors + sum(m4QAMdemod ~= message);
    % --------------------------------------------------- %
        end
end

% Check if we should modulate with 8PSK %
if (mux == [1, 1, 0, 0])
    errors = 0;
    m = 8;
    Es = 1;
    message = randi([0, m-1], [1, len_message]);
    for counter = 1:N
    % % 8PSK Modulation
        % According to item 1 - modulate, add noise, demodulate, accumulate error %
        m8PSK = pskmod(message, m);           
        noise = (sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) ... 
            + 1j*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        m8PSKnoisy = m8PSK + noise;
        m8PSKdemod = pskdemod(m8PSKnoisy, m); 
        errors = errors + sum(m8PSKdemod ~= message);
    % --------------------------------------------------- %
    end
end

% Check if we should modulate with 16-QAM %
if (mux == [1, 1, 1, 0])
    errors = 0;
    m = 16;
    Es = 5*2;
    message = randi([0, m-1], [1, len_message]);
    for counter = 1:N

    % % 16-QAM Modulation % %
        % According to item 1 - modulate, add noise, demodulate, accumulate error %
        m16QAM = qammod(message, m);            
        noise = (1)*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) ...
            + 1j*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        m16QAMnoisy = m16QAM + noise;
        m16QAMdemod = qamdemod(m16QAMnoisy, m); 
        errors = errors + sum(m16QAMdemod ~= message);
    % --------------------------------------------------- %
    end
end

% Check if we should modulate with 64-QAM %
if (mux == [1, 1, 1, 1])
    errors = 0;
    m = 64;
    Es = 42;
    message = randi([0, m-1], [1, len_message]);
    for counter = 1:N
    % % 64-QAM Modulation
        % According to item 1 - modulate, add noise, demodulate, accumulate error %
        m64QAM = qammod(message, m);
        noise = (1)*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message)) ...
            + 1j*(sqrt(Es/(2*(10^(EbN0vector(q)/10))*log2(m))))*randn(size(message));
        m64QAMnoisy = m64QAM + noise;
        m64QAMdemod = qamdemod(m64QAMnoisy, m); 
        errors = errors + sum(m64QAMdemod ~= message);
    % --------------------------------------------------- %
    end
end

% Add the found error count for this SNR %
Pe(q) = errors/(N*length(message));
end

%% Plotting Theoretical Curves %%
figure();
semilogy(EbN0vector(1: length(Pe)), Pe, 'black-o');       % Experimental Error Probability
hold on;

% Theoretical Modulation Probabilities %
semilogy(EbN0vector, Pe_tBPSK);
semilogy(EbN0vector, Pe_tQPSK);
semilogy(EbN0vector, Pe_t8PSK);
semilogy(EbN0vector, Pe_t16);
semilogy(EbN0vector, Pe_t64);

% Limiting Probability Line %
line([-2, 22], [probLim, probLim], 'LineStyle', '--', 'Color', 'black');

%%%%%%%% Band Assembly - Enable for Band Section %%%%%%%%%
% bands = [7.359999999999999, 10.960000000000001, 11.680000000000000, 16.240000000000002];
% line([bands(1), bands(1)], [10^-8, 1],  'LineStyle', '--', 'Color', 'black');
% line([bands(2), bands(2)], [10^-8, 1],  'LineStyle', '--', 'Color', 'black');
% line([bands(3), bands(3)], [10^-8, 1],  'LineStyle', '--', 'Color', 'black');
% line([bands(4), bands(4)], [10^-8, 1],  'LineStyle', '--', 'Color', 'black');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid on;
xlabel('Eb/No (dB)');
ylabel('SER');
title('Decided Modulation Curve');
axis([EbN0min, EbN0max, 10^-8, 1]);
legend('Experimental','BPSK','QPSK','8-PSK','16-QAM', '64-QAM');
hold off;