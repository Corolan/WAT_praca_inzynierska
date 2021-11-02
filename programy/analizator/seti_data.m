function seti_data()



tic; %rozpoczęcie stopera
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
warning('off');
disp("Hello!")
clear variables;%wyczyszczenie zmiennych przestrzeni roboczej

%deklaracje zmiennych globalnych
global fft_done;        % liczba wykonanych transformat
global Fs;              % Sampling frequency
global Ts;              % Sampling period
global center_freq;     %częstotliwość środkowa
global nmb_of_samples;  %liczba próbek
global analyze_threshold; %próg wykrywania pików
global gauss_analyze_threshold;%próg wykrywania krzywych gaussa
global tripulse_analyze_threshold;%próg wykrywania potrójnych impulsów
global nmb_of_pulses_threshold;
global bsmooth_boxcar_length;
global bsmooth_chunk_size;
global resol_thresh;

fft_done = 0;
%PARAMETRY SYGNAŁU
Fs= 9765.625;
Ts=1/Fs;
center_freq=1420019531.25;
%<spike_thresh>24</spike_thresh>
analyze_threshold = 24; 
%<gauss_peak_power_thresh>3.25</gauss_peak_power_thresh>
gauss_analyze_threshold = 3.25; % z pliku work_unit.sah
%gauss_analyze_threshold = 6.25; % z pliku work_unit.sah
%<triplet_thresh>9.73841</triplet_thresh>
tripulse_analyze_threshold = 9.738;
nmb_of_pulses_threshold=3;
if (nmb_of_pulses_threshold < 3)
    nmb_of_pulses_threshold=3;
end
nmb_of_samples=1048576;
bsmooth_boxcar_length=8192;% z pliku work_unit.sah
%pomnożone przez 32 by wygładzić całość w jednym przebiegu
bsmooth_chunk_size=32768;% z pliku work_unit.sah
resol_thresh=17000;% 17000
fake_switch=0;
pwr = [8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072]; %trans. fouriera będzie wykonywana dla fragmentów liczących tyle próbek

%pobranie i wstępne przygotowanie danych do analizy
disp('Wywołuję funkcję pobierania danych...');
data_raw = aquire_data;
disp('Dane pobrane');

disp('Wywołuje funkcję wygładzania danych...');
[data_after_smoothing, mean_value]=baseline_smoothing(data_raw, length(data_raw), bsmooth_boxcar_length, bsmooth_chunk_size, Fs, analyze_threshold, fake_switch);
disp('Wygładzanie wykonane');

%obliczenie transformat fouriera i zapisanie ich do plików
disp('Wywołuje funkcję obliczającą transformatę fouriera...');
fourier_transform(pwr, nmb_of_samples, data_after_smoothing, fft_done);
%fourier_transform(pwr, nmb_of_samples, data_raw, fft_done);
disp('Transformaty obliczone');

disp('Wywołuję funkcję rysującą wykres 3D...');
%funkcja rysuje wykresy 3D
%draw_surface(pwr, Ts, Fs, analyze_threshold, nmb_of_samples);
%disp('Wykresy narysowane');

%analiza danych
disp('Wywołuję funkcję analizy danych...');

disp('Wywołuję funkcję wykrywania pików...');
peak_searching(pwr, analyze_threshold, Fs, mean_value);
disp('Wykryte piki zapisane do pliku');

disp('Wykrywanie charakterystyk gaussowskich...');
gauss_searching(pwr, gauss_analyze_threshold, resol_thresh, Fs, Ts, nmb_of_samples, mean_value);
disp('Wykryte dopasowania zapisane do pliku');

disp('Wykrywanie regularnych impulsów...');
tripulse_searching(pwr, nmb_of_pulses_threshold, tripulse_analyze_threshold, Ts, Fs, nmb_of_samples, mean_value);
disp('Wykryte ciągi trzech impulsów zapisane do pliku');

disp('Koniec funkcji analizy danych!');
disp('***Program zakończył pracę***');

toc; %koniec stopera
end

% definicje funkcji
function z = aquire_data()%wczytanie danych z pliku
filename = sprintf('wynik-5.txt');
delimiterIn = ' ';
raw_data_from_file = importdata(filename,delimiterIn);
clear filename delimiterIn;
Re = raw_data_from_file(:,1);
Im = raw_data_from_file(:,2);
z=complex(Re, Im)';% utworzenie liczb zespolonych z pobranych kolumn danych TRANSPOZYCJA
clear raw_data_from_file Re Im;%wyczyszczenie niepotrzebnych zmiennych
end
function [data_after_smoothing, mean_value]=baseline_smoothing(z, data_length, bsmooth_boxcar_length, bsmooth_chunk_size, Fs, analyze_threshold, fake_switch)
% %funkcja dzieli dane na kawałki o długości bsmooth_chunk_size, następnie
% %wygładza je za pomocą średniej ruchomej...
% num_of_smooths = data_length/bsmooth_chunk_size;
% frequency_domain_before_smoothing = fft(z,data_length);
% for cnt=1:1:num_of_smooths
%     part_before_smoothing=frequency_domain_before_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt); %aktualnie przetwarzany kawałek danych
%     tmp_bs=(smooth(part_before_smoothing, bsmooth_boxcar_length, 'moving')); 
%     %... i składa je z powrotem w jeden wektor
%     frequency_domain_after_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt)=tmp_bs;
% end
    %inny sposób - za  pomocą fukcji filter
    num_of_smooths = data_length/bsmooth_chunk_size;
    frequency_domain_before_smoothing = fft(z,data_length);
    f=linspace(0, Fs, data_length);
    xl='częstotliwość [Hz]';
    yl='Moduł';
    
    fake_sig_info=' ';
    if (fake_switch == 1)
        fake_sig_info=' z fałszywym sygnałem';
        %*******FAKE_SIGNAL****************
        szer_piku=400;
        signal_pwr=2150; %2000!!!
        %250;
        offset=263000;
        %460000;
        %463000;
        %fake_signal_abs=[zeros(1,(data_length-szer_piku)/2-offset), (signal_pwr)*ones(1,szer_piku), zeros(1,(data_length-szer_piku)/2+offset)];
        fake_signal_abs=[zeros(1,offset), (signal_pwr)*ones(1,szer_piku), zeros(1,(data_length-szer_piku-offset))];
        
        %t=linspace(0,1,data_length);
        %fake_signal_abs_gauss=signal_pwr*sin(1.10*t);
        
        %fake_signal=(fake_signal_abs);
        fake_signal=hilbert(fake_signal_abs);
        fake_signal_abs_value=(abs(fake_signal));
        disp('Hello');
        fprintf('\nMax, wartosc sygnału sztucznego: %d\n',max(fake_signal_abs_value));
        
        disp('Hello-v2');
                
        figure; c=plot(f,abs(fake_signal));
        title('Fałszywy sygnał')
        xlabel(xl);
        ylabel(yl);
        saveas(c, 'fakesignal', 'jpeg')
        
        %fake_signal_time=ifft(fake_signal, length(fake_signal));
        [~, I]=max(abs(fake_signal_abs).^2);
        %max_val=sprintf('\nMax: %d Index: %d\n',M,I)
%         max_val_freq=sprintf('\nCzęstotliwość sygnału sztucznego: %d Hz',f(I));
%         disp(max_val_freq);
        fprintf('\nCzęstotliwość sygnału sztucznego: %d Hz\n',f(I));
    figure;
    before_smooth_plot=plot(f, abs(frequency_domain_before_smoothing));
    plot_name=sprintf('Dane przed wygładzeniem\n(dziedzina częstotliwości) \n');
    xlabel(xl);
    ylabel(yl);
    title(plot_name);
    saveas(before_smooth_plot, 'Dane_przed_wygładzeniem-freq-no-fake', 'jpeg');
    comp1=abs(frequency_domain_before_smoothing); %porównanie fake-nofake
        %****************
        %frequency_domain_after_smoothing=frequency_domain_after_smoothing+fake_signal;
%         figure;
%         before_smooth_fake_signal=plot(f, abs(frequency_domain_before_smoothing).^2);
%         plot_name=sprintf('Dane po wygładzeniu + sztuczny sygnał\n o częstotliwości: %d Hz\n(dziedzina częstotliwości)', f(I));
%         hold on;
%         plot(f, mean_value*ones(1, data_length), 'r','LineWidth',2);
%         plot(f, analyze_threshold*mean_value*ones(1, data_length), 'g','LineWidth',2);
%         legend('Dane','Wartość średnia', 'Próg detekcji');
%         xlabel(xl);
%         ylabel(yl);
%         title(plot_name);
%         hold off;
%         saveas(after_smooth_fake_signal, plot_name, 'jpeg');
        frequency_domain_before_smoothing=frequency_domain_before_smoothing+fake_signal;
        comp2=abs(frequency_domain_before_smoothing); %porównanie fake-nofake
    end
    
    

    figure;
    before_smooth_plot=plot(f, abs(frequency_domain_before_smoothing));
    plot_name=sprintf('Dane przed wygładzeniem\n(dziedzina częstotliwości) \n%s', fake_sig_info);
    xlabel(xl);
    ylabel(yl);
    title(plot_name);
    saveas(before_smooth_plot, 'Dane_przed_wygładzeniem-freq', 'jpeg');
    
    %parametry filtra
    windowSize = bsmooth_boxcar_length; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;

    frequency_domain_after_smoothing=(filter(b,a,frequency_domain_before_smoothing));
    
       
    
%     %%%FILTER w kawałkach%%%
%     for cnt=1:1:num_of_smooths
%         part_before_smoothing=frequency_domain_before_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt); %aktualnie przetwarzany kawałek danych
%         tmp_bs=filter(b, a, part_before_smoothing);
%         %... i składa je z powrotem w jeden wektor
%         frequency_domain_after_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt)=tmp_bs;
%     end
    mean_value = (sum(abs(frequency_domain_after_smoothing).^2))/length(frequency_domain_after_smoothing);
    %odzyskanie danych w dziedzinie czasu
    
    
    figure;
    after_smooth_plot=plot(f, abs(frequency_domain_after_smoothing).^2);
    plot_name=sprintf('Dane po wygładzeniu\n(dziedzina częstotliwości)');
%     hold on;
%     plot(f, mean_value*ones(1, data_length), 'r','LineWidth',2);
%     plot(f, analyze_threshold*mean_value*ones(1, data_length), 'g','LineWidth',2);
%     legend('Dane','Wartość średnia', 'Próg detekcji');
    xlabel(xl);
    ylabel('Kwadrat modułu');
    title(plot_name);
    hold off;
    saveas(after_smooth_plot, plot_name, 'jpeg');
    data_after_smoothing=ifft(frequency_domain_after_smoothing, data_length);
    
    if (fake_switch == 1)
        fprintf('\nWzgledna wartosc sygnału sztucznego: %d razy poziom średni\n',max(fake_signal_abs_value)/mean_value);
        comp=abs(comp1-comp2);
        figure; cmp_plot=plot(f, comp);
        xlabel(xl);
        ylabel(yl);
        axis ([-inf inf -inf 12000]);
        title('Porównanie sygnału przed i po wygładzeniu');
        saveas(cmp_plot, 'porownanie', 'jpeg');
    end
    
    
    

end
function fourier_transform(pwr, nmb_of_samples, data_after_smoothing, fft_done)
% w poniższej pętli program oblicza transformaty F. zadanej długości zapisując je do plików tekstowych

%all_ffts=fft_all_count(pwr,nmb_of_samples);%liczba transformat do wykonania
for k=1:1:length(pwr) %kolejne kroki obliczeń
    
    %fprintf('Powinno być 2^17, czyli 131072: %d\n', pwr(k));
    nmb_of_fft=nmb_of_samples/(pwr(k)); %liczba transformat niezbędnych do wykonania w danym kroku
    act_data_length = nmb_of_samples/nmb_of_fft;%liczba próbek w aktualnie badanym fragmencie danych
    
    %prealokacja pamięci dla macierzy zawierającej wykonane fft
    fft_matrix_current_step=zeros(nmb_of_fft, act_data_length);
    %fprintf('Zalokowana pamięć: %d wierszy %d kolumn\n', nmb_of_fft, act_data_length);
    
    for cnt=1:nmb_of_fft % cnt - numer obecnie obliczanej transformaty
        
        %wycinanie kawałków danych
        act_data=data_after_smoothing((nmb_of_samples / nmb_of_fft)*(cnt-1)+1:(nmb_of_samples / nmb_of_fft)*cnt); %aktualnie przetwarzany kawałek danych
        %obliczanie transformaty fouriera aktualnego wycinka danych
        tmp_fft=fft(act_data, length(act_data));%transformata fouriera badanego fragmentu
        fft_value = (abs(tmp_fft)).^2;
        %norm_factor=max(abs(tmp_fft));
        %fft_value = abs(tmp_fft)/norm_factor;%normalizacja
        fft_done = fft_done+1;
        %wpisanie obliczonej transformaty do macierzy
        fft_matrix_current_step(cnt, 1:act_data_length) = fft_value;
    end
    current_fft_matrix_label= sprintf('fft_matrix/length_%d.txt', pwr(k));
    
    norm_factor = max(max(fft_matrix_current_step));%unormowanie z całej macierzy
    fft_matrix_current_step=fft_matrix_current_step/norm_factor;
    dlmwrite(current_fft_matrix_label, fft_matrix_current_step, 'newline', 'pc');
    
end
disp('Wykonano wszystkie transformaty i zapisano rezultaty do plików');
fprintf('Wykonano następującą liczbę transformat: %d\n', fft_done)

end
function peak_searching(pwr, analyze_threshold, Fs, mean_value)
plots_drawn=0;

file_peak_result = fopen('file_peak_result.txt','w');
for k=1:1:length(pwr)
    
    %wczytanie pliku z danymi
    actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
    actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
    [nmb_of_fft, size_of_fft]=size(actual_fft_matrix);
    for j=1:nmb_of_fft %wybieram pojedyncze wiersze z macierzy - pojedyncze fft
        actual_fft_row=actual_fft_matrix(j, :);
        
        %poprzedni sposob wyszukiwania
        %actual_mean_value=mean(actual_fft_row);
        
        %nowy sposób: suma wartości próbek (energia) dzielona przez liczbę próbek 
        %actual_mean_value=sum(actual_fft_row)/length(actual_fft_row);
    
        %values_indices_over_threshold=find(actual_fft_row >= (analyze_threshold*actual_mean_value));
        values_indices_over_threshold=find(actual_fft_row >= (analyze_threshold*mean_value));
        
        if(isempty(values_indices_over_threshold) == 0)%jeżeli macierz wyników powyzej threshold jest nie pusta
            f=linspace(0,Fs,pwr(k));
            %zapisanie do pliku informacji o wykrytym piku:
            for p_count=1:length(values_indices_over_threshold)
                fprintf(file_peak_result,'Plik: length_%d.txt, nr wiersza: %d, częstotliwość: %d, Wartość unormowana: %d/1\n',pwr(k), j, f(values_indices_over_threshold(p_count)), actual_fft_matrix(j, values_indices_over_threshold(p_count)));
            end
            figure;%!!!!!
            fig_act_fft = stairs(f, actual_fft_row, 'b');
            hold on;
            plot(f, mean_value*ones(1,size_of_fft), 'r'); % dorysowanie wartości średniej do wykresu
            plot(f, analyze_threshold*mean_value*ones(1,size_of_fft), 'r'); % dorysowanie wartości granicznej do wykresu
            %linie wartości średniej i granicznej
            %mean_value_label=sprintf('Wartość średnia: %03.2f',actual_mean_value);
            %text(f(1), 1.2*actual_mean_value, mean_value_label,'FontSize',7);
            %thresh_value_label=sprintf('Próg: %03.2f Wartość graniczna: %03.2f', analyze_threshold, analyze_threshold*actual_mean_value);
            %text(f(1), 0.9*analyze_threshold*actual_mean_value, thresh_value_label,'FontSize',7);
            
            plot_title=sprintf('Widmo sygnału dla transformaty nr %d o długości %d próbek\n', j, size_of_fft);
            xyz=sprintf('Wartość średnia: %d Wartość progowa: %d', mean_value,  analyze_threshold*mean_value);
            %plot_subtitle=text(xyz);
            %plot_subtitle.FontSize =7;
            title( {plot_title; xyz},'FontSize',7);
            %axis([-inf inf 0 1.1]);
 
            set(gca,'XTickLabelRotation',45)
            xlabel('f (Hz)');
            ylabel('Widmo amplitudowe');
            
            %wyłączenie wypisywania wartości przekraczających próg na
            %wykres
%             for n=1:length(values_indices_over_threshold)
%                 max_value_label=sprintf('%2.2f', actual_fft_matrix(j, values_indices_over_threshold(n)));
%                 text(f(values_indices_over_threshold(n)), actual_fft_matrix(j, values_indices_over_threshold(n)), max_value_label,'FontSize',5);
%             end
            hold off;
            plots_drawn=plots_drawn+1;
            %fig_file_name= sprintf('plots/plots_%d/fig_%d_%d', pwr(k), size_of_fft, j);
            fig_file_name= sprintf('test/%d-fig_%d_%d', plots_drawn, size_of_fft, j);
            saveas(fig_act_fft, fig_file_name, 'jpeg');
            
        end
    end
end
fprintf('narysowano %d wykres(ów)\n', plots_drawn);
fclose(file_peak_result);
end
function gauss_searching(pwr, gauss_analyze_threshold, resol_thresh, Fs, Ts, nmb_of_samples, mean_value)

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
        %kolumny mają po 12 sekund długości
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
                
                %actual_mean_value=sum(actual_freq_col)/length(actual_freq_col);
                %actual_mean_value=mean(actual_freq_col);
                
                values_indices_over_threshold=find(actual_freq_col>=(gauss_analyze_threshold*mean_value));
                [h, ~, stats] = chi2gof(actual_freq_col, 'cdf', @(z)normcdf(z, mean(actual_freq_col)));
                gauss_searching_done=gauss_searching_done+1; %
                %if(isempty(values_indices_over_threshold) == 0 && h == 0 && stats.chi2stat > 3.5)%jeżeli macierz wyników powyzej threshold jest nie pusta

                if(isempty(values_indices_over_threshold) == 0 && h == 0 && stats.chi2stat > 4.0)%jeżeli macierz wyników powyzej threshold jest nie pusta
                    if (sum(values_indices_over_threshold)>=length(actual_freq_col)/4 && sum(values_indices_over_threshold)<=(3*length(actual_freq_col))/4)
                    
                        freq=Fs*col_cnt/pwr(k);
                        fprintf(file_gauss_result,'Transf. o dłg: %d , fragment: %d/%d, nr kolumny: %d, częstotliwość: %d\n', pwr(k), d, D_chunk, col_cnt, freq);
                        t=linspace(0,D_Z*pwr(k)*Ts,length(actual_freq_col));
                        fig_act_gauss = plot(t,actual_freq_col, 'b');
                        fig_act_gauss.LineWidth = 1.0;
                        
                        hold on;
                        %plot_title=sprintf('Dopasowanie gaussa dla częstotliwości\n %d dla transformaty o długości %d punktów %d\nchi2stat= %d', col_cnt, pwr(k), nmb_of_fft, stats.chi2stat);
                        actual_thresh_value=gauss_analyze_threshold*actual_mean_value;
                        plot_title=sprintf('Dopasowanie gaussa nr %d/%d/%d \n częstotliwość %d [Hz]\nWartość średnia: %d Wartość progu: %d\n', d, col_cnt, pwr(k), freq, actual_mean_value, actual_thresh_value);
                        title(plot_title);
                        
                        mline=plot(t, actual_mean_value*ones(1,length(actual_freq_col)), 'r-.'); % dorysowanie wartości średniej do wykresu
                        mline.LineWidth = 1.2;
                        %mean_label=sprintf('%d', actual_mean_value);
                        %text(0, 0.95*actual_mean_value, mean_label, 'Color','black');
                        tline=plot(t, gauss_analyze_threshold*actual_mean_value*ones(1,length(actual_freq_col)), 'g:'); % dorysowanie wartości granicznej do wykresu
                        tline.LineWidth = 1.2;
                        %thresh_label=sprintf('%d', gauss_analyze_threshold*actual_mean_value);
                        %text(0, 0.95*gauss_analyze_threshold*actual_mean_value, thresh_label, 'Color','black');
                        fig_file_name= sprintf('gauss_fit/%d_fig_%d_%d_%d',gauss_plots_done, d, col_cnt, pwr(k));
                        f = fit(t.',actual_freq_col.' ,'gauss1');
                        p=plot(f, 'black');
                        p.LineWidth = 1.5;
                        legend('Dane', 'Wartość średnia', 'Wartość progu', 'Dopasowana krzywa', 'Location', 'northeast','Orientation','vertical');
                        ylabel('Widmo amplitudowe');
                        xlabel('t [s]');
                        hold off;
                        saveas(fig_act_gauss, fig_file_name, 'jpeg');
                        gauss_plots_done=gauss_plots_done+1;
                    end
                end
            end
            
            
        end
    end
end
fprintf('Wykonano %d poszukiwań gaussa\n', gauss_searching_done);
fprintf('Narysowano %d dopasowań gaussa\n', gauss_plots_done);
fclose(file_gauss_result);
end
function draw_surface(pwr, Ts, Fs, analyze_threshold, nmb_of_samples)
    for k=1:1:length(pwr)
        actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
        actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
        
        fig_act_fft=surf(actual_fft_matrix,'EdgeColor','none');
        [row_nmb, col_nmb] = size(actual_fft_matrix);
        title_label=sprintf('Wykres macierzy transformaty\n o długości %d punktów', col_nmb);
        title(title_label);
        fig_act_fft.XData=linspace(0,Fs,col_nmb);
        fig_act_fft.YData=linspace(0,Ts*nmb_of_samples,row_nmb);
        %view([0, 90]);
        
        xlabel('f [Hz]');
        ylabel('t [s]');
        set(gca,'YDir','reverse');
        zlabel('widmo ampl.');
        fig_file_name=sprintf('fft_matrix/surface_%d', pwr(k));
        %saveas(fig_act_fft,fig_file_name, 'fig');
        saveas(fig_act_fft,fig_file_name, 'png');
        view([0, 90]);%view([7, 70]);
        fig_file_name=sprintf('fft_matrix/surface_%d_cv', pwr(k));
        saveas(fig_act_fft,fig_file_name, 'png');
    end
% for k=1:1:length(pwr)
%     actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
%     actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
%     [nmb_of_fft, size_of_fft]=size(actual_fft_matrix);%[liczba wierszy, liczba kolumn]
%     for j=1:nmb_of_fft %wybieram pojedyncze kolumny z macierzy - pojedyncze częstotliwości
%         actual_freq_row=(actual_fft_matrix(j, :));
%         actual_mean_value=mean(actual_freq_row);
%     end
%     
% end
end
function total_nmb_of_tripulse=tripulse_counter(pwr)
tmp=0;
for k=1:1:length(pwr)
    %wczytanie pliku z danymi
    actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
    actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
    [~, size_of_fft]=size(actual_fft_matrix);
    tmp=tmp+size_of_fft;
    total_nmb_of_tripulse=tmp;
end
end
function tripulse_searching(pwr, nmb_of_pulses_threshold, tripulse_analyze_threshold, Ts, Fs, nmb_of_samples, mean_value)
file_tripulse_result = fopen('file_tripulse_result.txt','w');
total_nmb_of_tripulse=tripulse_counter(pwr);%zliczenie wszystkich fitów do zrobienia
fprintf('Do sprawdzenia %d częstotliwości na obecność pulsów\n', total_nmb_of_tripulse);
tripulse_searching_done=0;
tripulse_plots_done=0;
for k=1:1:length(pwr)
    %wczytanie pliku z danymi
    actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
    actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
    [nmb_of_fft, size_of_fft]=size(actual_fft_matrix);%[liczba wierszy, liczba kolumn]

    for j=1:size_of_fft %wybieram pojedyncze kolumny z macierzy - pojedyncze częstotliwości
        actual_freq_col=(actual_fft_matrix(:, j))';

        %actual_mean_value=mean(actual_freq_col);
        %actual_mean_value=sum(actual_freq_col)/length(actual_freq_col);
        
        values_indices_over_threshold=find(actual_freq_col >= (tripulse_analyze_threshold*mean_value));
        if(length(values_indices_over_threshold) >= nmb_of_pulses_threshold)%jeżeli macierz wyników powyzej threshold zawiera co najmniej trzy elementy
            %poniższy fragment kodu sprawdza, czy wartości powyżej progu leżą w równych odległościach od siebie
            interval=values_indices_over_threshold(2)-values_indices_over_threshold(1);
            for h=1:(length(values_indices_over_threshold)-1)
                tmp_interval=values_indices_over_threshold(h+1)-values_indices_over_threshold(h);
                if(tmp_interval==interval)
                    continue; %jeśli interwał jest stały, przejdź do następnej wartości
                else
                    interval=0;
                    break; % jeśli jest zmienny, nie szukaj dalej
                end
            end
            
            if(interval ~= 0) %jeśli interwał pozostał stały dla wszystkich wartości z danego wektora czestotliwości
                %rysuj wykres
                t = linspace(0,Ts*nmb_of_samples,length(actual_freq_col)); %wektor czasu
                time_resol=Ts*nmb_of_samples/length(actual_freq_col);
                fig_act_tripulse = plot(t,actual_freq_col, 'b');
                for g=1:length(values_indices_over_threshold)
                    text(t(values_indices_over_threshold(g)), actual_freq_col(values_indices_over_threshold(g)), '@','FontSize',7);
                end
                axis([-inf inf -inf inf]);
                hold on;
                plot(t, actual_mean_value*ones(1,nmb_of_fft), 'r'); % dorysowanie wartości średniej do wykresu
                plot(t, tripulse_analyze_threshold*actual_mean_value*ones(1,nmb_of_fft), 'r'); % dorysowanie wartości granicznej do wykresu
                %linie wartości średniej i granicznej
                mean_value_label=sprintf('Wartość średnia: %03.2f',actual_mean_value);
                %text(t(1), actual_mean_value, mean_value_label,'FontSize',7);
                thresh_value_label=sprintf('Próg: %03.2f   Wartość graniczna: %03.2f', tripulse_analyze_threshold, tripulse_analyze_threshold*actual_mean_value);
                legend('Dane', mean_value_label, thresh_value_label);
                %text(t(1), tripulse_analyze_threshold*actual_mean_value, thresh_value_label,'FontSize',7);
                freq=Fs*j/pwr(k);
                plot_title=sprintf('Potrójne impulsy dla częstotliwości %d [Hz]\ntransformaty o długości %d punktów\nInterwał: %d [s]', freq, pwr(k), interval*time_resol);
                title(plot_title);
                set(gca,'XTickLabelRotation',45)
                xlabel('t [s]');
                ylabel_label=sprintf('Wartości transformowanych próbek\nw kolejnych chwilach czasu');
                ylabel(ylabel_label);
                hold off;
                fig_file_name= sprintf('tripulse_search/fig_%d_%d', nmb_of_fft, j);
                saveas(fig_act_tripulse, fig_file_name, 'jpeg');
                tripulse_plots_done=tripulse_plots_done+1;
                %zapisanie do pliku informacji o znalezionych impulsach
                fprintf(file_tripulse_result,'Plik: length_%d.txt, nr kolumny: %d, częstotliwość: %d\n',pwr(k), j, freq);
            end
        end
        tripulse_searching_done=tripulse_searching_done+1;
    end
end
fprintf('Wykonano %d poszukiwań impulsów\n', tripulse_searching_done);
fprintf('Narysowano %d poszukiwań impulsów\n', tripulse_plots_done);
fclose(file_tripulse_result);
end


%{
function gauss_searching_backup(pwr, gauss_analyze_threshold, resol_thresh, Fs, Ts, nmb_of_samples)
gauss_searching_done=0;
gauss_plots_done=0;
file_gauss_result = fopen('file_gauss_result.txt','w');
for k=1:1:length(pwr)
    if pwr(k)<=resol_thresh %ograniczenie rozdzielczości
        %wczytanie pliku z danymi
        actual_fft_matrix_label=sprintf('fft_matrix/length_%d.txt', pwr(k));
        actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
        [nmb_of_fft, size_of_fft]=size(actual_fft_matrix);%[liczba wierszy, liczba kolumn]
        
        %podzielenie aktualnej macierzy danych na macierz, której wiersze
        %mają po 12 sekund długości
        D_Z=floor(12/pwr(k)*Fs); %liczba transformat o długości pwr(k) zapewniająca 12 sekund sygnału
        D_chunk=(nmb_of_samples/pwr(k)-mod(nmb_of_samples/pwr(k),D_Z))/D_Z;
        
        actual_fft_matrix( (1+nmb_of_fft-mod(nmb_of_samples/pwr(k),D_Z)):nmb_of_fft, :)=[];
        %mam już przyciętą macież, teraz trzeba podzielić ją na D_chunk
        %kawałków
        
        %dzielenie na D_chunk kawałków
        for d=1:1:D_chunk
            actual_gauss_matrix=actual_fft_matrix( 1+(d-1)*D_Z : d*D_Z, :);
            %ANALIZA teraz należy przeanalizować poszczególne kolumny
            %macierzy actual_gauss_matrix zawiera ona 12 sekundowe
            %fragmenty sygnału
            for col_cnt=1:size_of_fft
                actual_freq_col=(actual_gauss_matrix(:, col_cnt)');
                actual_mean_value=mean(actual_freq_col);
                values_indices_over_threshold=find(actual_freq_col>=(gauss_analyze_threshold*actual_mean_value));
                [h, ~, stats] = chi2gof(actual_freq_col, 'cdf', @(z)normcdf(z, mean(actual_freq_col)));
                if(isempty(values_indices_over_threshold) == 0 && h == 0 && stats.chi2stat > 3)%jeżeli macierz wyników powyzej threshold jest nie pusta
                    if (sum(values_indices_over_threshold)>length(actual_freq_col)/2)
                        freq=Fs*col_cnt/pwr(k);
                        fprintf(file_gauss_result,'Transf. o dłg: %d , fragment: %d/%d, nr kolumny: %d, częstotliwość: %d\n', pwr(k), d, D_chunk, col_cnt, freq);
                        t=linspace(0,D_Z*pwr(k)*Ts,length(actual_freq_col));
                        fig_act_gauss = plot(t,actual_freq_col, 'b');
                        
                        hold on;
                        %plot_title=sprintf('Dopasowanie gaussa dla
                        %częstotliwości\n %d dla transformaty o długości %d
                        %punktów %d\nchi2stat= %d', col_cnt, pwr(k),
                        %nmb_of_fft, stats.chi2stat);
                        plot_title=sprintf('Dopasowanie gaussa nr %d/%d/%d \n częstotliwość %d [Hz]', d, col_cnt, pwr(k), freq);
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
end


function fourier_transform_backup(pwr, nmb_of_samples, data_after_smoothing, fft_done)
% w poniższej pętli program oblicza transformaty F. zadanej długości zapisując je do plików tekstowych

%all_ffts=fft_all_count(pwr,nmb_of_samples);%liczba transformat do wykonania
for k=1:1:length(pwr) %kolejne kroki obliczeń
    
    fprintf('Powinno być 2^17, czyli 131072: %d\n', pwr(k));
    nmb_of_fft=nmb_of_samples/(pwr(k)); %liczba transformat niezbędnych do wykonania w danym kroku
    act_data_length = nmb_of_samples/nmb_of_fft;%liczba próbek w aktualnie badanym fragmencie danych
    
    %prealokacja pamięci dla macierzy zawierającej wykonane fft
    fft_matrix_current_step=zeros(nmb_of_fft, act_data_length);
    fprintf('Zalokowana pamięć: %d wierszy %d kolumn\n', nmb_of_fft, act_data_length);
    
    for cnt=1:nmb_of_fft % cnt - numer obecnie obliczanej transformaty
        
        %wycinanie kawałków danych
        act_data=data_after_smoothing((nmb_of_samples / nmb_of_fft)*(cnt-1)+1:(nmb_of_samples / nmb_of_fft)*cnt); %aktualnie przetwarzany kawałek danych
        %obliczanie transformaty fouriera aktualnego wycinka danych
        tmp_fft=fft(act_data);%transformata fouriera badanego fragmentu
        fft_value = (abs(tmp_fft));
        %norm_factor=max(abs(tmp_fft));
        %fft_value = abs(tmp_fft)/norm_factor;%normalizacja
        fft_done = fft_done+1;
        %wpisanie obliczonej transformaty do macierzy
        fft_matrix_current_step(cnt, 1:act_data_length) = fft_value;
    end
    current_fft_matrix_label= sprintf('fft_matrix/length_%d.txt', pwr(k));
    
    norm_factor = max(max(fft_matrix_current_step));
    fft_matrix_current_step=fft_matrix_current_step/norm_factor;
    dlmwrite(current_fft_matrix_label, fft_matrix_current_step, 'newline', 'pc');
    
end
disp('Wykonano wszystkie transformaty i zapisano rezultaty do plików');
fprintf('Wykonano następującą liczbę transformat: %d\n', fft_done)

end

function data_after_smoothing_backup=baseline_smoothing(z, data_length, bsmooth_boxcar_length, bsmooth_chunk_size)
% %funkcja dzieli dane na kawałki o długości bsmooth_chunk_size, następnie
% %wygładza je za pomocą średniej ruchomej...
% num_of_smooths = data_length/bsmooth_chunk_size;
% frequency_domain_before_smoothing = fft(z,data_length);
% for cnt=1:1:num_of_smooths
%     part_before_smoothing=frequency_domain_before_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt); %aktualnie przetwarzany kawałek danych
%     tmp_bs=(smooth(part_before_smoothing, bsmooth_boxcar_length, 'moving')); 
%     %... i składa je z powrotem w jeden wektor
%     frequency_domain_after_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt)=tmp_bs;
% end
%odzyskanie danych w dziedzinie czasu

%inny sposób - za  pomocą fukcji filter
frequency_domain_before_smoothing = fft(z,data_length);
windowSize = bsmooth_boxcar_length; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
filtered_data=(filter(b,a,frequency_domain_before_smoothing));
data_after_smoothing=ifft(filtered_data, data_length);
end

%}
