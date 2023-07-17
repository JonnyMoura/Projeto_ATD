%ex1
close all;
clear all;

labels = load("labels.txt");

colors = ["red" "blue" "green" "cyan" "magenta" "yellow" "#55FF55" "#F506FF" "#7E2F8E" "#4DBEEE" "#D95319" "#ABC120"]; % cor associada respetivamente a cada atividade

activ_label = ["WALKING", "WALKING_UPSTAIRS", "WALKING_DOWNSTAIRS", "SITTING", "STANDING", "LAYING", "STAND_TO_SIT", "SIT_TO_STAND", "SIT_TO_LIE", "LIE_TO_SIT", "STAND_TO_LIE", "LIE_TO_STAND"]; % atividades

experiment = ["acc_exp42_user21.txt", "acc_exp43_user21.txt", "acc_exp44_user22.txt", "acc_exp45_user22.txt", "acc_exp46_user23.txt","acc_exp47_user23.txt", "acc_exp48_user24.txt", "acc_exp49_user24.txt"]; % ficheiros dos acelerometros

eixos = ["x" "y" "z"];
%ex 2
for k=1:length(experiment)
    user = k + 41;

    indices = find(labels(:, 1) == user); % vai guardar o indice do label desse user a cada iteraçao
    linhas = labels(indices, 3:end); % vai guardar a atividade e os valores das 2 ultimas colunas do label file

    file = load(experiment(k));

    Fs = 50;
    Ts = 1/Fs;
    N = length(file);
    t = 0: Ts: N*Ts-Ts;

    figure() % para o eixo dos x
    subplot(3,1,1);
    plot(t, file(:,1), "black-");
    hold on;

    minValX = min(file(:,1));
    maxValX = max(file(:,1));
    
    subplot(3,1,2); % para o eixo dos y
    plot(t, file(:,2), "black-");
    hold on;
    
    minValY = min(file(:,2));
    maxValY = max(file(:,2));

    subplot(3,1,3); % para o eixo dos z
    plot(t, file(:,3), "black-");
    hold on;

    minValZ = min(file(:,3));
    maxValZ = max(file(:,3));
    
    for a=1:length(linhas)
        maximo = linhas(a, 3);

        minimo = linhas(a, 2);

        atividade = linhas(a, 1);

        tt = t(minimo:maximo);

        subplot(3,1,1);
        p = plot(t(minimo:maximo), file(minimo:maximo,1));
        p.Color = colors(atividade);
        
        ac_tag = activ_label(atividade);
        ac_tag_X = tt(1);
        ac_tag_Y = minValX + 0.1;
        if (mod(a, 2) == 0)
            ac_tag_Y = 0.9*maxValX;
        end 
        text(ac_tag_X, ac_tag_Y, ac_tag, 'FontSize', 5);
        
        hold on;

        subplot(3,1,2);
        p = plot(t(minimo:maximo), file(minimo:maximo,2));
        p.Color = colors(atividade);

        ac_tag = activ_label(atividade);
        ac_tag_X = tt(1);
        ac_tag_Y = minValY + 0.1;
        if (mod(a, 2) == 0)
            ac_tag_Y = 0.9*maxValY;
        end  
        text(ac_tag_X, ac_tag_Y, ac_tag, 'FontSize', 5);
        
        hold on;

        subplot(3,1,3);
        p = plot(t(minimo:maximo), file(minimo:maximo,3));
        p.Color = colors(atividade);

        ac_tag = activ_label(atividade);
        ac_tag_X = tt(1);
        ac_tag_Y = minValZ + 0.1;
        if (mod(a, 2) == 0)
            ac_tag_Y = 0.9*maxValZ;
        end 
        text(ac_tag_X, ac_tag_Y, ac_tag, 'FontSize', 5);
        
        hold on;
    
    end    
end

%3.1/3.2

num_exp = [42 43 44 45 46 47 48 49]; % numero das experiencias

freqs = cell(1,12); % frequencias relevantes
ExpTable = cell(1, 8); % experiencias
amps = cell(1,12); % amplitudes relevantes

ExpTable{8} = load(experiment(8));
ExpTable{7} = load(experiment(7));
ExpTable{6} = load(experiment(6));
ExpTable{5} = load(experiment(5));
ExpTable{4} = load(experiment(4));
ExpTable{3} = load(experiment(3));
ExpTable{2} = load(experiment(2));
ExpTable{1} = load(experiment(1));

for a = 1:8 % percorrer as experiencias
    for k = 1:12 % para cada atividade
            lbsExp = labels(labels(:,1) == num_exp(a), :);
            [amps{a,k}, freqs{a,k}] = calculo_dft(ExpTable{a}, lbsExp, k, 50); % guarda nas cells as amplitudes e frequencias (relevantes)
    end
end

%ex 3.3

for a = 1:3 % percorre as 3 atividades dinamicas walking, ...
    for e = 1:3 % isto para cada eixo (x, y, z)
        P_tot = 0;
        for k = 1:numel(freqs{1,a}{e})
            P_tot = freqs{1,a}{e}(k) + P_tot;
        end

        fprintf("desvio padrao: %s -> %s", activ_label(a), eixos(e));
        disp(std([freqs{1,a}{:,e}]));
        
        fprintf("media: %s -> %s", activ_label(a), eixos(e));
        disp(P_tot*60/numel(freqs{1,a}{:,e}));
        
    end
end

%ex 3.4
figure();
hold on
for i = 1:8
    for k = 1:11
        legend(activ_label(1), activ_label(2), activ_label(3), activ_label(4), activ_label(5), activ_label(6),activ_label(7), activ_label(8), activ_label(9), activ_label(10), activ_label(11), activ_label(12));
        simbolo = ['+', '+', '+', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*'];
        plot3(freqs{i,k}{1}(1), freqs{i,k}{2}(1), freqs{i,k}{3}(1), simbolo(k));
        view(3);
    end
end
hold off

%ex 3.5
simbolo = ['*', '-', 'o', '*', '-', 'o', '*', '-', 'o', 'p', 's', 'x'];

% para as atividades dinamicas
figure();
hold on
for i = 1:8
    for k = 1:3
        plot3(amps{i,k}{1}(1), amps{i,k}{2}(1), amps{i,k}{3}(1), simbolo(k));
        legend(activ_label(1), activ_label(2), activ_label(3));
        view(3);
    end
end
hold off

% para as atividades estaticas
figure();
hold on
for i = 1:8
    for k = 4:6
        plot3(amps{i,k}{1}(1), amps{i,k}{2}(1), amps{i,k}{3}(1), simbolo(k));
        legend(activ_label(4), activ_label(5), activ_label(6));
        view(3);
    end
end
hold off

% para as atividades de transiçao
figure();
hold on
for i = 1:8
    for k = 7:11
        plot3(amps{i,k}{1}(1), amps{i,k}{2}(1), amps{i,k}{3}(1), simbolo(k));
        legend(activ_label(7), activ_label(8), activ_label(9), activ_label(10), activ_label(11), activ_label(12));
        view(3);
    end
end
hold off

%ex 4.1

d = load(experiment(1)); % usar apenas a experiencia 42

data = d(:,3); % ir buscar apenas o eixo dos z da experiencia

expLabels = labels(labels(:,1) == 42, :);

janelas_dft(d, expLabels, activ_label);

%ex 4.2 + 4.3
plotSTFT(data, 50)

function [amp, freq] = calculo_dft(data, labels, atividade, Fs)

    amp = {};
    freq = {};
    
    exp = find(labels(:,3) == atividade);
    
    n = numel(exp);
    
    for i = 1:n
        
        m = exp(i);
        fim = labels(m,5);
        inicio = labels(m, 4);
        
        intervalo = fim - inicio + 1;
        deltfreq = Fs/intervalo;
        
        if mod(intervalo, 2) ~= 0
            f = -Fs/2+deltfreq/2:deltfreq:Fs/2-deltfreq/2;
        else
            f = -Fs/2:deltfreq:Fs/2-deltfreq;
         
        end
        
        janela = inicio:fim;
        eixo_x = data(janela, 1);
        eixo_y = data(janela, 2);
        eixo_z = data(janela, 3);
        
        eixo_x = detrend(eixo_x);    
        eixo_y = detrend(eixo_y);    
        eixo_z = detrend(eixo_z);
            
        DFT_x = abs(fftshift(fft(eixo_x))); % calcular o dft do eixo dos x
        DFT_y = abs(fftshift(fft(eixo_y))); % calcular o dft do eixo dos y
        DFT_z = abs(fftshift(fft(eixo_z))); % calcular o dft do eixo dos z
        
        magnitude_x = abs(DFT_x); % calcular a magnitude do eixo dos x
        magnitude_y = abs(DFT_y); % calcular a magnitude do eixo dos y
        magnitude_z = abs(DFT_z); % calcular a magnitude do eixo dos z
        
        threshold_x = 0.8*max(magnitude_x);
        threshold_y = 0.8*max(magnitude_y);
        threshold_z = 0.8*max(magnitude_z);
        
        [peaks_x,locs_x] = findpeaks(magnitude_x, 'MinPeakHeight',threshold_x);
        freq_x = f(locs_x);
        freq_x = freq_x(freq_x > 0);
        
        [peaks_y,locs_y] = findpeaks(magnitude_y, 'MinPeakHeight',threshold_y);
        freq_y = f(locs_y);
        freq_y = freq_y(freq_y > 0);
        
        [peaks_z,locs_z] = findpeaks(magnitude_z, 'MinPeakHeight',threshold_z);
        freq_z = f(locs_z);
        freq_z = freq_z(freq_z > 0);
        
        amp{1} = peaks_x;
        amp{2} = peaks_y;
        amp{3} = peaks_z;

        freq{1} = freq_x;
        freq{2} = freq_y;
        freq{3} = freq_z;
        
    end
end 

function [] = plotSTFT(data, Fs)

N = numel(data);
Ts = 1/Fs;
t = N*Ts;

intervalos_janela = [0.2 0.1 0.05 0.001 0.007];

for i = 1:5
    T_part = intervalos_janela(i)*t;
    
    T_Overlap = T_part/2;
    
    N_part = round(T_part*Fs); % frame
    
    janela = blackman(N_part); % usou-se a janela de blackman
    
    Noverlap = round(T_Overlap*Fs);
    
    f = linspace(-Fs/2,Fs/2,N_part);
    x =  find(f>=0);

    espetro = [];
    
    for j = 1:N_part-Noverlap:N-N_part
        x_part = data(j:j+N_part-1).*janela; % frame
        magn_x_part = abs(fftshift(fft(x_part)));
        espetro = horzcat(espetro,magn_x_part(x));
    end
    figure();
    imagesc(20*log10(espetro))
    
    set(gca,'YDir','normal')
    xlabel('t(s)')
    ylabel('f(Hz)')
    title(['STFT | tam janela: ',num2str(intervalos_janela(i), '%.3fxt | OverLap de 50%')])
    colorbar
end
end

function janelas_dft(dados, labels, atividades, dyn, show_wv, wv_I)
    arguments
        dados
        labels
        atividades 
        dyn = true
        show_wv = false
        wv_I = 1
    end

    counts = zeros(1,12);
    Fs = 50;
    
    for j = 1:size(labels, 1)
        num_ativ = labels(j,3);
        if (num_ativ ~= 1)
            continue
        end
        ativ = string(atividades(num_ativ));
        inicio = labels(j, 4);
        fim = labels(j,5);
        intervalo = fim - inicio + 1;
        delt_freq = Fs/intervalo;
        
        if mod(intervalo, 2) == 0
            f = -Fs/2:delt_freq:Fs/2-delt_freq;
        else
            f = -Fs/2+delt_freq/2:delt_freq:Fs/2-delt_freq/2;
        end
        
        figure('Name', sprintf("%s - Janelas p/ DFT", ativ));
        
        % janela Retangular
        janela = inicio:fim;
        eixo_z = dados(janela, 3);
        dft_eixo_z = abs(fftshift(fft(eixo_z)));
        
        subplot(4,1,1)
        plot(f, dft_eixo_z)
        title("dft eixo z")
        xlabel("freq(Hz)")
        ylabel("magn")
        
        % janela de Hamming 
        janela_hamming = hamming(intervalo);
        z_janelado = eixo_z.*janela_hamming;
        
        dft_eixo_z = abs(fftshift(fft(z_janelado)));
        
        subplot(4,1,2)
        plot(f, dft_eixo_z)
        title("dft eixo z janela Hamming")
        xlabel("freq(Hz)")
        ylabel("magn")
        
        % janela de Hann
        janela_hann = hann(intervalo);
        z_janelado = eixo_z.*janela_hann;
        dft_eixo_z = abs(fftshift(fft(z_janelado)));
        
        subplot(4,1,3)
        plot(f, dft_eixo_z)
        title("dft eixo z janela Hann")
        xlabel("freq(Hz)")
        ylabel("magn")
        
        % janela de Blackman
        janela_blackman = blackman(intervalo);
        z_janelado = eixo_z.*janela_blackman;
        
        dft_eixo_z = abs(fftshift(fft(z_janelado)));
        
        subplot(4,1,4)
        plot(f, dft_eixo_z)
        title("dft eixo z janela Blackman")
        xlabel("freq(Hz)")
        ylabel("magn")
        
        counts(num_ativ) = 1;
        if (show_wv == true) && (j == wv_I)
            wvtool(janela_hamming, janela_hann, janela_blackman);
        end
        return
    end
end
