%% Exame Final de EET-50: Luis Felipe Silva Rezende Soares %%
% Instituto Tecnológico de Aeronáutica, Junho de 2023 %
% Introdução a Sistemas de Comunicação %

%% Questão 1 %%

% Estime, atraves de simulacao computacional (m ́etodo de Monte Carlo),
% a probabilidade de erro de sımbolo (Pe) em fun ̧c ̃ao da raz ̃ao Eb/N0, para
% um canal AWGN de media zero e PSD N0/2, para as modulacoes em banda
% base a seguir:
% a) 2-PAM (BPSK)
% b) 4-PAM
% c) 4-QAM (QPSK)
% d) 2-FSK
% Obs: Os resultados das simula ̧c ̃oes devem ser comparados com os
% resultados te ́oricos, na mesma figura de Pe × (Eb/N0)

%% Resolução da Questão 1 %%
clear all; clc;

% Calculo da energia de simbolo%
delta = 1;
M = 2;                          % Modulacao 2-PAM
N = 4;                          % Modulacao 4-PAM
S = 4;                          % Modulacao 4-QAM
F = 2;                          % Modulacao 2-FSK
Es = (delta^2)*(M^2 - 1)/12;    % Energia do Simbolo de 2-PAM
Es_4PAM = (delta^2)*(N^2-1)/12; % Energia do Simbolo de 4-PAM
Es_4QAM = (delta^2)*(S-1)/6;    % Energia do Símbolo de 4-QAM
Es_2FSK = (delta^2);            % Energia do Símbolo de 2-FSK (='s)

% Calculo da energia de bit %

K_2PAM = ceil(log(M)/log(2));           % Ceiling de numero fracionario
Eb0 = Es/K_2PAM;                        % Energia do Bit de 2PAM
K_4PAM = ceil(log(N)/log(2));           % Ceiling de numero fracionario de 4PAM
Eb0_4PAM = Es_4PAM/K_4PAM;              % Energia do Bit de 4PAM
K_4QAM = ceil(log2(S));                 % Ceiling do Bit de 4QAM
Eb0_4QAM = Es_4QAM/K_4QAM;              % Energia do Bit de 4QAM
K_2FSK = ceil(log2(F));                 % Ceiling do Bit de 2FSK
Eb0_2FSK = Es_2FSK/K_2FSK;              % Energia do Bit de 2FSK


% Parametros Iniciais de Calculo $
mu = 0;
sigma = 1;
delta_n = 10^-2;

% Inicializacao dos Vetores e Valores %
Q_v = [];                       % Vetor para as integrais em Eb/N0;
Q_v_4PAM = [];                  % Vetor para as integrais em Eb/N0 de 4-PAM;
Q_v_4QAM = [];
Q_v_2FSK = [];
P_erro_2PAM_v = [];             % Vetor para as probabilidades de erro 2PAM
P_erro_4PAM_v = [];             % Vetor para as probabilidades de erro 4PAM
P_erro_4QAM_v = [];
P_erro_2FSK_v = [];
Ebv = [];                       % Vetor para os valores de Energia de Bit de 2PAM
Ebv_4PAM = [];                  % Vetor para os valores de Energia de Bit de 4PAM
Ebv_4QAM = [];
Ebv_2FSK = [];
P_e_2PAM_DB = [];               % Vetor para os valores de DB de 2PAM
P_e_4PAM_DB = [];               % Vetor para os valores de DB de 4PAM
P_e_4QAM_DB = [];
P_e_2FSK_DB = [];
N0 = 2*Eb0;                     % A razão Eb/N0 inicia em -3 dB, em que Eb é igual à PSD do ruido AWGN.
IterLimit = 100*Eb0/delta_n + 1; % Limite de Interações: 30 dB %

%% Loop de Construção dos Vetores %%

% Iteração para ambas as modulacoes %
for i = 1:IterLimit

    Eb = Eb0 + i*delta_n;                       % Atualizacao da Energia do Bit de 2-PAM
    Eb_4PAM = Eb0_4PAM + i*delta_n;             % Atualizacao da Energia do Bit de 4-PAM
    Eb_4QAM = Eb0_4QAM + i*delta_n;             % Atualizacao da Energia do Bit de 4-QAM
    Eb_2FSK = Eb0_2FSK + i*delta_n;             % Atualizacao da Energia do Bit de 2-FSK
    Ebv(i) = 10*log10(Eb/N0);                   % 2-PAM em DB
    Ebv_4PAM(i) = 10*log10(Eb_4PAM/N0);         % 4-PAM em DB    
    Ebv_4QAM(i) = 10*log10(Eb_4QAM/N0);         % 4-QAM em DB
    Ebv_2FSK(i) = 10*log10(Eb_2FSK/N0);         % 2-FSK em DB

% O calculo da probabilidade de erro de uma modulação M-PAM é dada por:
x = sqrt((6*(log2(M))/(M^2-1))*(Eb/N0));
lowerLimit_2PAM = x;
y = sqrt((6*(log2(N))/(N^2-1))*(Eb_4PAM/N0));
lowerLimit_4PAM = y;
upperLimit = Inf;

% Definindo a Funcao Gaussiana de Densidade de Probabilidade %
q = @(x) normpdf(x, mu, sigma);

% Calculando a integral da Gaussiana para Modulacao 2-PAM %
Q = integral(q, lowerLimit_2PAM, upperLimit);    % Integral
Q_v(i) = Q;                                 % Push em Vetor de 2-PAM
P_erro_2PAM = 2*(1 - (1/M))*Q;              % Probabilidade de 2-PAM
P_erro_2PAM_v(i)= P_erro_2PAM;              % Push no Vetor de Probabilidade
P_e_2PAM_DB(i) = 10*log10(P_erro_2PAM);     % Push no Vetor de Probabilidade em DB

% Calculando a integral da Gaussiana para Modulacao 4-PAM %
Q = integral(q, lowerLimit_4PAM, upperLimit);
Q_v_4PAM(i) = Q;
P_erro_4PAM = 2*(1 - (1/N))*Q;
P_erro_4PAM_v(i) = P_erro_4PAM;
P_e_4PAM_DB(i) = 10*log10(P_erro_4PAM);

%  % Vamos agora calcular para uma Modulacao 4-QAM %  %

% Calculando a integral da Gaussiana para Modulacao 4-QAM
lowerLimit_4QAM = sqrt((3*log2(S)/(S-1))*(Eb_4QAM/N0));
Q = integral(q, lowerLimit_4QAM, upperLimit);
Q_v_4QAM = Q;
P_erro_4QAM = 4*(1 - (1/sqrt(S)))*Q - 4*((1 - (1/sqrt(S)))^2)*Q^2;
P_erro_4QAM_v(i) = P_erro_4QAM;
P_e_4PAM_DB(i) = 10*log10(P_erro_4QAM);

% % Vamos calcular para a Modulacao 2-FSK % %
lowerLimit_2FSK = sqrt(Eb_2FSK/N0);
Q = integral(q, lowerLimit_2FSK, upperLimit);
Q_v_2FSK = Q;
P_erro_2FSK = Q;
P_erro_2FSK_v(i) = P_erro_2FSK;
P_e_2FSK_DB(i) = 10*log10(P_erro_2FSK);
end

%% Plot da Figura com Todas as Modulacoes %%
figure();
semilogy(Ebv, P_erro_2PAM_v);
hold on;
xlabel('Eb/N0 (dB)');
ylabel('SER');
title('SER vs Eb/N0');
semilogy(Ebv, P_erro_4PAM_v);
semilogy(Ebv, P_erro_4QAM_v)
semilogy(Ebv, P_erro_2FSK_v);
ylim([10^-6, 0.5]);
legend('2-PAM', '4-PAM','4-QAM','2-FSK', 'Location', 'southwest');
grid on;
hold off;