%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This Source file compares the Bit Error Rate for 4-QAM digital modulation of two different constelation mappings:    %
%  The mapping done via Natural Mapping (00, 01, 10, 11), done via vector [0, 3, 2, 1], and Gray (00, 01, 11, 10)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;
clear all;

% 4-QAM Modulation %
m = 4;
messageLength = 100;     % Message Length
N = 3000;               % Number of Monte Carlo Iterations
quantity = 50;

% Iteration Parameters % 
EbNomin = -2;
EbNomax = 16;
binaryErrorCounter = 0;
grayErrorTotal = [];
grayErrorCounter = 0;
probError4QAMbinary = [];
probError4QAMgray = [];
EbN0vector = EbNomin:(EbNomax - EbNomin)/quantity:EbNomax;

% Calculation of Symbol Energy: It is twice the sqrt(4) - PAM, thus twice
% the 2-PAM
Es = 2;

% Let's plot the behavior of the modulation as Eb/N0 increases.
for q = 1:quantity

    grayError = 0;
    binaryError = 0;

    for i = 1:N

        % Message Generation and its Binary Representation %
        message = randi([0, m-1], [1, messageLength]);                                  
        messageBin = dec2bin(message);
   
        SNR = EbN0vector(q);     
       
        % Modulation of the Message with Noise. % 
        messageQAMbinary = qammod(message, m, [0, 3, 2, 1]);    
        messageQAMgray = qammod(message, m, 'gray'); 

        % Adding Noise (in both real and complex parts) to the modulated signal %
        binaryNoise = (sqrt(Es/(2*(10^(SNR/10))*log2(m))))*randn(size(message)) ... 
            + 1j*(sqrt(Es/(2*(10^(SNR/10))*log2(m))))*randn(size(message));        
        grayNoise = (sqrt(Es/(2*(10^(SNR/10))*log2(m))))*randn(size(message)) ... 
            + 1j*(sqrt(Es/(2*(10^(SNR/10))*log2(m))))*randn(size(message)); 
        messageQAMnoisyBinary = messageQAMbinary + binaryNoise; 
        messageQAMnoisyGray = messageQAMgray + grayNoise;    
        
        % QAM Demodulation %
        messageDemodQAMbinary = qamdemod(messageQAMnoisyBinary, m, [0, 3, 2, 1]);
        messageDemodQAMgray = qamdemod(messageQAMnoisyGray, m, 'gray');
        messageDemodQAMgrayBinary = dec2bin(messageDemodQAMgray);
        messageDemodQAMbinaryBinary = dec2bin(messageDemodQAMbinary);
    
        % Error Computation (both in Gray and Natural Mapping) %
        grayError = grayError + sum(sum(messageDemodQAMgrayBinary ~= messageBin));
        binaryError = binaryError + sum(sum(messageDemodQAMbinaryBinary ~= messageBin));
    end

    % Adding to the respective probability vectors %
    probError4QAMbinary(q) = (binaryError/(2*N))/(messageLength);
    probError4QAMgray(q) = (grayError/(2*N))/(messageLength);
end

%% 4-QAM Curve Plot %%

m = 4;                                                                      

SNRlinear = 10.^(EbN0vector./10);

% Calculate the theoretical error probability
argument = sqrt(3*(log2(m)/(m-1)).*(SNRlinear));
 
% Probability Function of 4-QAM Modulation %
 Pe_t = (1/2)*4*(1 - 1/sqrt(m)).*qfunc(argument) - 4*(1 - 1/sqrt(m))^2.*qfunc(argument).^2; % Q Function

% Theoretical curve plotted%
figure();
semilogy(EbN0vector, Pe_t,'black');
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Error Probability');
title('4-QAM Natural vs 4-QAM Gray Error');
semilogy(EbN0vector(1:length(probError4QAMbinary)), probError4QAMbinary);
semilogy(EbN0vector(1:length(probError4QAMgray)), probError4QAMgray);
axis([-2, 10, 10^-5, 1]);
legend('Gray Theoretical Approximation', 'Natural Mapping', 'Gray Mapping');
hold off;