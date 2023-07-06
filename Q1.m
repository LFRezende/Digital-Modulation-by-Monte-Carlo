%% Exame Final de EET-50: Luis Felipe Silva Rezende Soares %%
% Instituto Tecnológico de Aeronáutica, Junho de 2023 %
% Introdução a Sistemas de Comunicação

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utilized for the Final Exam of Communication Systems (EET-50)) in the Aeronautics Institute of Technology,ITA, Brazil %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ### ITEM A ### --- Obtencao do Probabilidade de Erro %%
format long;

% Mensagem %
clear all;
clc;
close all;

% Parametros de Iteracao % 
EbN0min = -2;
EbN0max = 16;
quantidade = 100;
errorCounter = 0;
probErro2Pam = [];
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantidade:EbN0max;

% Modulacao 2PAM:
m = 2;
Es = 1;                 % Energia de Simbolo Es = (1^2 + 1^2)/2 = 1;
N = 1000;               % Iteracao em Monte Carlo
len_mensagem = 1280;    % Comprimento da Mensagem enviada (Simbolos)
mensagem = randi([0, m-1], [1, len_mensagem]);  % Mensagem

% Vamos plotar o comportamento da modulacao conforme aumenta-se o Eb/N0.%

for q=1:quantidade
    % Contador de Erros Totais Obtidos apos Monte Carlo %
    errorCounter = 0; 
    %%% Laco da Rotina Iterativa de Monte Carlo %%%
    for i= 1:N
        
        % Selecao de Eb/N0 para ruido construido %
        EbN0 = EbN0vector(q);                               
        mensagemPAM = pammod(mensagem, m);     % Modulacao PAM  

        % Construcao do Ruido por Energia de Simbolo e Adicao %
        ruido = (sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(mensagem));
        mensagemPAMruidosa = mensagemPAM + ruido;

        % Demodulacao PAM %
        mensagemDemodPAM = pamdemod(mensagemPAMruidosa, m);             
    
        % Analise: Diferencas entre o sinal recebido e o sinal original.%
        errorCounter = errorCounter + sum(mensagemDemodPAM ~= mensagem);
    end
    % Computo de Probabilidade%
    probErro2Pam(q) = errorCounter/(len_mensagem*N);
end

% Plot da Teorica 2-PAM %%

% Parametros Basicos
m = 2; % m-PAM
% Conjunto de Dados para o Plot da Teorica % 

EbN0linear = 10.^(EbN0vector./10);

% Probabilidade teorica 
argumento = sqrt(6*(log2(m)/(m^2-1)).*(EbN0linear));
Pe_t = 2*(1 - 1/m)*qfunc(argumento); % Formula Teorica de Probabilidade

% Plot da Figura
figure();
semilogy(EbN0vector, Pe_t, 'b');
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Probabilidade de Erro');
title('Probabilidade de Erro 2-PAM');
semilogy(EbN0vector(1:length(probErro2Pam)), probErro2Pam, 'r');
axis([-2, 10, 10^-5, 1]);
legend('Teórico', 'Experimental');
hold off;



%% ### ITEM B ### --- Calculo da Modulacao 4-PAM %%

% Parametros de Iteracao % 
EbN0min = -2;
EbN0max = 16;
quantidade = 100;
errorCounter = 0;
probErro2Pam = [];
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantidade:EbN0max;
probErro4Pam = [];
% Modulacao 4PAM:
m = 4;
mensagem = randi([0, m-1], [1, len_mensagem]);
% Energia de Símbolo: 3^2 + 1^2 + 1^2 + 3^2 = 20 --> Es = 20/4.
Es = 5;  % Energia do Simbolo

% Vamos plotar o comportamento da modulacao conforme aumenta-se o Eb/N0.
for q = 1:quantidade
    errorCounter = 0;
    for i= 1:N
                                                                       
        EbN0 = EbN0vector(q);                                                              
        
        % Normalização da Mensagem e feita no ruido. %
        mensagemPAM = pammod(mensagem, m);
        ruido = (sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(mensagem));
        mensagemPAMruidosa = mensagemPAM + ruido;
        mensagemDemodPAM = pamdemod(mensagemPAMruidosa, m);
    
        % Analise: Diferencas entre o sinal recebido e o sinal original.%
        errorCounter = errorCounter + sum(mensagemDemodPAM ~= mensagem);     
    end
    probErro4Pam(q) = errorCounter/(len_mensagem*N);
end



% Plot da Teorica 4-PAM %

% Parametros Basicos
m = 4; % m-PAM
EbN0linear = 10.^(EbN0vector./10);

% Calcular a probabilidade de erro teorica
argumento = sqrt(6*(log2(m)/(m^2-1)).*(EbN0linear));

% Normalizacao pela Energia de Simbolo.
%Eb = Es/K => eS = KEb --> M-pam: 
Pe_t = 2*(1 - 1/m)*qfunc(argumento); % Funcao Q

% Curva teorica plotada
figure();
semilogy(EbN0vector, Pe_t, 'black');
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Probabilidade de Erro');
title('Probabilidade de Erro da 4-PAM');
semilogy(EbN0vector(1:length(probErro4Pam)), probErro4Pam, 'green');
axis([-2, 14, 10^-5, 1]);
legend('Teórico', 'Experimental');
hold off;

%% ### ITEM C ###  --- Calculo da Probabilidade de Erro de Modulacao QAM %%


% Parametros de Iteracao % 
EbN0min = -2;
EbN0max = 16;
errorCounter = 0;
probErro4QAM = [];
% Modulacao 4-QAM:
m = 4;
mensagem = randi([0, m-1], [1, len_mensagem]);
% Cálculo da Energia de Símbolo: É o dobro da sqrt(4) - PAM, logo o dobro
% da 2-PAM
Es = 2;

% Vamos plotar o comportamento da modulacao conforme aumenta-se o Eb/N0.
for q=1:quantidade
    errorCounter = 0;  
    for i= 1:N
        EbN0 = EbN0vector(q); 

        % Modulacao, ruido e demodulacao %
        mensagemQAM = qammod(mensagem, m); 
        ruido = (sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(mensagem)) ... 
            + 1j*(sqrt(Es/(2*(10^(EbN0/10))*log2(m))))*randn(size(mensagem));
        mensagemQAMruidosa = mensagemQAM + ruido; 
        mensagemDemodQAM = qamdemod(mensagemQAMruidosa, m);
    
        % Analise: Diferencas entre o sinal recebido e o sinal original.%
        errorCounter = errorCounter + sum(mensagemDemodQAM ~= mensagem);

    end
    probErro4QAM(q) = errorCounter/(len_mensagem*N);
end


% Plot da Curva 4-QAM %

                                                                   
EbN0linear = 10.^(EbN0vector./10);

% Calcular a probabilidade de erro teorica
argumento = sqrt(3*(log2(m)/(m-1)).*(EbN0linear));
Pe_t = 4*(1 - 1/sqrt(m)).*qfunc(argumento) - 4*(1 - 1/sqrt(m))^2.*qfunc(argumento).^2; % Funcao Q

% Curva teorica plotada
figure();
semilogy(EbN0vector, Pe_t);
hold on;
grid on;
xlabel('Eb/No (dB)');
ylabel('Probabilidade de Erro');
title('Probabilidade de Erro da 4-QAM');
semilogy(EbN0vector(1:length(probErro4QAM)), probErro4QAM);
axis([-2, 10, 10^-5, 1]);
legend('Teórico', 'Experimental');
hold off;


%% ### ITEM D ###  - Modulacao FSK %%
m = 2;         % Modulacao 2-FSK
k = log2(m);   % Bits
EbNo = 5;      % Eb/No (dB)
Fs = 16;       % Taxa de Amostragem
nsamp = 8;     % Numero de Amostras por Simbolo
freqsep = 10;  % Separacao em Frequencia 

len_mensagem = 128;
N = 100;
mensagem = randi([0, m-1], [1, len_mensagem]);


% Parametros de Iteracao % 
EbN0min = -2;
EbN0max = 16;
probErro2FSK = [];
EbN0vector = EbN0min:(EbN0max - EbN0min)/quantidade:EbN0max;
Es = 1; % A energia de simbolo do sinal BFSK é a propria amplitude do simbolo.

% Modulacao 2FSK
for q=1:quantidade
    taxaDeErro = 0;
    for i= 1:N
    
        EbN0 = EbN0vector(q);                                                                 
        mFSK = fskmod(mensagem,m,freqsep,nsamp,Fs);
        mFSKRuidosa  = awgn(mFSK,(EbN0)+10*log10(k)-10*log10(nsamp),'measured',[],'dB');
        mFSKDemod = fskdemod(mFSKRuidosa,m,freqsep,nsamp,Fs);

        %   Para a modulacao 2FSK, a existencia de 1 bit torna BER = SER %
        [num,BER] = biterr(mensagem,mFSKDemod);

        taxaDeErro = taxaDeErro + BER;
        
        
    end
    %   Probabilidade de Erro Experimental %
    probErro2FSK(q) = taxaDeErro/N;

    %   Vetor de Probabilidade de Erro Teorico (em 2-FSK, BER = SER)%
    BER_teorico= berawgn(EbN0,'fsk',m,'noncoherent');
    Pe_t(q) = BER_teorico;
end

% Plot da Curva para 2-FSK
figure();
semilogy(EbN0vector(1:length(Pe_t)), Pe_t);
hold on;
semilogy(EbN0vector(1:length(probErro2FSK)), probErro2FSK);
grid on;
xlabel('Eb/No (dB)');
ylabel('Probabilidade de Erro');
title('Probabilidade de Erro da 2-FSK');
legend('Teórico', 'Experimental');
hold off;
