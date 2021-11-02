tic; %rozpoczęcie stopera
clc;
%set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
warning('off');
disp("Hello!");
clear variables;%wyczyszczenie zmiennych przestrzeni roboczej

%deklaracje zmiennych globalnych
global fft_done;    % liczba wykonanych transformat
global Fs;            % Sampling frequency                    
global Ts;             % Sampling period
global center_freq; %częstotliwość środkowa
global nmb_of_samples;
global gauss_analyze_threshold;

%PARAMETRY SYGNAŁU
fft_done = 0;
Fs= 9765.625;% z pliku work_unit.sah
Ts=1/Fs;
center_freq=1420019531.25;% z pliku work_unit.sah
%<gauss_peak_power_thresh>3.25</gauss_peak_power_thresh>
gauss_analyze_threshold = 3.25; % z pliku work_unit.sah
nmb_of_samples=1048576;


pwr = [8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072]; %trans. fouriera będzie wykonywana dla fragmentów liczących tyle próbek

resol_thresh=17000;% 17000
    
    gauss_searching_done=0;
    gauss_plots_done=0;
file_gauss_result = fopen('file_gauss_result.txt','w');
    for k=1:1:length(pwr)
        if pwr(k)<=resol_thresh %ograniczenie rozdzielczości
            %wczytanie pliku z danymi
            actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
            actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
            [nmb_of_fft, size_of_fft]=size(actual_fft_matrix);%[liczba wierszy, liczba kolumn]
            
            %podzielenie aktualnej macierzy danych na macierz, której
            %wiersze mają po 12 sekund długości
            D_Z=floor(12/pwr(k)*Fs); %liczba transformat o długości pwr(k) zapewniająca 12 sekund sygnału
            D_chunk=(nmb_of_samples/pwr(k)-mod(nmb_of_samples/pwr(k),D_Z))/D_Z;
            
            actual_fft_matrix( (1+nmb_of_fft-mod(nmb_of_samples/pwr(k),D_Z)):nmb_of_fft, :)=[];
            %mam już przyciętą macież, teraz trzeba podzielić ją na D_chunk kawałków
            
            %dzielenie na D_chunk kawałków
            for d=1:1:D_chunk
                actual_gauss_matrix=actual_fft_matrix( 1+(d-1)*D_Z : d*D_Z, :);
                %ANALIZA
                %teraz należy przeanalizować poszczególne kolumny macierzy actual_gauss_matrix
                %zawiera ona 12 sekundowe fragmenty sygnału
                for col_cnt=1:size_of_fft
                    actual_freq_col=(actual_gauss_matrix(:, col_cnt)');
                    actual_mean_value=mean(actual_freq_col);
                    values_indices_over_threshold=find(actual_freq_col>=(gauss_analyze_threshold*actual_mean_value));
                    [h, p, stats] = chi2gof(actual_freq_col, 'cdf', @(z)normcdf(z, mean(actual_freq_col)));
                    if(isempty(values_indices_over_threshold) == 0 && h == 0 && stats.chi2stat > 3)%jeżeli macierz wyników powyzej threshold jest nie pusta
                        fprintf(file_gauss_result,'Transf. o dłg: %d , fragment: %d/%d, nr kolumny: %d', pwr(k), d, D_chunk, col_cnt);
                        if (sum(values_indices_over_threshold)>length(actual_freq_col)/2)
                            t=linspace(0,D_Z*pwr(k)*Ts,length(actual_freq_col));
                            fig_act_gauss = plot(t,actual_freq_col, 'b');

                            hold on;
                            %plot_title=sprintf('Dopasowanie gaussa dla częstotliwości\n %d dla transformaty o długości %d punktów %d\nchi2stat= %d', col_cnt, pwr(k), nmb_of_fft, stats.chi2stat);
                            plot_title=sprintf('Dopasowanie gaussa nr %d/%d/%d ', d, col_cnt, pwr(k));
                            title(plot_title);

                            plot(t, actual_mean_value*ones(1,length(actual_freq_col)), 'r'); % dorysowanie wartości średniej do wykresu
                            mean_label=sprintf('%d', actual_mean_value);
                            text(0,0.95*actual_mean_value, mean_label, 'Color','red');
                            plot(t, gauss_analyze_threshold*actual_mean_value*ones(1,length(actual_freq_col)), 'g'); % dorysowanie wartości granicznej do wykresu
                            thresh_label=sprintf('%d', gauss_analyze_threshold*actual_mean_value);
                            text(0,0.97*gauss_analyze_threshold*actual_mean_value, thresh_label, 'Color','red');
                            fig_file_name= sprintf('gauss_fit/fig_%d_%d_%d', d, col_cnt, pwr(k));
                            f = fit(t.',actual_freq_col.','gauss1');
                            plot(f, 'black');
                            legend('Dane', 'Wartość średnia', 'Wartość progu', 'Dopasowana krzywa', 'Location', 'bestoutside','Orientation','horizontal');
                            ylabel('Widmo amplitudowe');
                            xlabel('t [s]');
                            hold off;
                            saveas(fig_act_gauss, fig_file_name, 'jpeg');
                            gauss_plots_done=gauss_plots_done+1;
                        end
                    end
                end
                
            gauss_searching_done=gauss_searching_done+1;    
            end
        end
    end
     fprintf('Wykonano %d poszukiwań gaussa\n', gauss_searching_done);
     fprintf('Narysowano %d dopasowań gaussa\n', gauss_plots_done);
fclose(file_gauss_result);
toc;
