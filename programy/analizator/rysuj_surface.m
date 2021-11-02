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

for g=1:5
    fprintf("Przebieg nr %d", g);
    %pobranie i wstępne przygotowanie danych do analizy
    disp('Wywołuję funkcję pobierania danych...');
    data_raw = aquire_data(g);
    disp('Dane pobrane');
    disp('Wywołuje funkcję wygładzania danych...');
    [data_after_smoothing, ~]=baseline_smoothing(data_raw, length(data_raw), bsmooth_boxcar_length, bsmooth_chunk_size, Fs, analyze_threshold, fake_switch, g);
    disp('Wygładzanie wykonane');
    %obliczenie transformat fouriera i zapisanie ich do plików
    disp('Wywołuje funkcję obliczającą transformatę fouriera...');
    %fourier_transform(pwr, nmb_of_samples, data_after_smoothing, fft_done, g);
    disp('Transformaty obliczone');
    disp('Wywołuję funkcję rysującą wykres 3D...');
    %draw_surface(pwr, Ts, Fs, analyze_threshold, nmb_of_samples, g);
    disp('Wykresy narysowane');
end


toc; %koniec stopera
end

% definicje funkcji
function z = aquire_data(g)%wczytanie danych z pliku
filename = sprintf('wynik-%d.txt',g);
delimiterIn = ' ';
raw_data_from_file = importdata(filename,delimiterIn);
clear filename delimiterIn;
Re = raw_data_from_file(:,1);
Im = raw_data_from_file(:,2);
z=complex(Re, Im)';% utworzenie liczb zespolonych z pobranych kolumn danych TRANSPOZYCJA
clear raw_data_from_file Re Im;%wyczyszczenie niepotrzebnych zmiennych
end
function [data_after_smoothing, mean_value]=baseline_smoothing(z, data_length, bsmooth_boxcar_length, bsmooth_chunk_size, Fs, analyze_threshold, fake_switch, g)
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
    
    %fake_sig_info=' ';   
    

%     figure;
%     before_smooth_plot=plot(f, abs(frequency_domain_before_smoothing)/max(abs(frequency_domain_before_smoothing)));
%     plot_name=sprintf('Dane przed wygładzeniem\n(dziedzina częstotliwości)-plik nr %d', g);
%     file_plot_name=sprintf('Dane_przed_wygładzeniem_plik_%d', g);
%     xlabel(xl);
%     ylabel(yl);
%     title(plot_name);
%     saveas(before_smooth_plot, file_plot_name, 'jpeg');
    

    
    %parametry filtra
    windowSize = bsmooth_boxcar_length; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;

    frequency_domain_after_smoothing=(filter(b,a,frequency_domain_before_smoothing));
    
       
    
%     %%%FILTER w kawałkach%%%
%     for cnt=1:1:num_of_smooths
%         part_before_smoothing=frequency_domain_before_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt); %aktualnie przetwarzany kawałek danych
%         %tmp_bs=filter(b, a, part_before_smoothing);
%         tmp_bs=smooth(part_before_smoothing,bsmooth_boxcar_length);
%         %... i składa je z powrotem w jeden wektor
%         frequency_domain_after_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt)=tmp_bs;
%     end
    
    
    %frequency_domain_after_smoothing=(filter(b,a,frequency_domain_before_smoothing));
    y=abs(frequency_domain_after_smoothing).^2;
    mean_value = (sum(abs(frequency_domain_after_smoothing).^2))/length(frequency_domain_after_smoothing)/max(y);
    %odzyskanie danych w dziedzinie czasu
    
    
%     figure;
%     
%     after_smooth_plot=plot(f, y/max(y) );
%     plot_name=sprintf('Dane po wygładzeniu\n(dziedzina częstotliwości) - plik %d',g);
%     file_plot_name=sprintf('Dane_po_wygładzeniu_plik_%d',g);
%     hold on;
%     plot(f, mean_value*ones(1, data_length), 'r','LineWidth',2);
%     plot(f, analyze_threshold*mean_value*ones(1, data_length), 'g','LineWidth',2);
%     legend('Dane','Wartość średnia', 'Próg detekcji');
%     xlabel(xl);
%     ylabel('Kwadrat modułu');
%     title(plot_name);
%     hold off;
%     saveas(after_smooth_plot, file_plot_name, 'jpeg');


     data_after_smoothing=ifft(frequency_domain_after_smoothing, data_length);
    
    figure;
    t=linspace(0, data_length/Fs, data_length);
    after_smooth_plot=plot(t, abs(data_after_smoothing));
    plot_name=sprintf('Dane po wygładzeniu\n(dziedzina czasu)-plik nr %d', g);
    file_plot_name=sprintf('Dane_po_wygładzeniem_czas_plik_%d', g);
    axis ([0, data_length/Fs, -inf, inf]);
    xlabel("czas [s]");
    ylabel("");
    title(plot_name);
    saveas(after_smooth_plot, file_plot_name, 'jpeg');
         

end
function fourier_transform(pwr, nmb_of_samples, data_after_smoothing, fft_done, g)
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
        norm_factor=max(fft_value);
        fft_value = fft_value/norm_factor;%normalizacja
        fft_done = fft_done+1;
        %wpisanie obliczonej transformaty do macierzy
        fft_matrix_current_step(cnt, 1:act_data_length) = fft_value;
    end
    current_fft_matrix_label= sprintf('fft_matrix/test/length_%d_plik_%d.txt', pwr(k),g);
    
    %norm_factor = max(max(fft_matrix_current_step));%unormowanie z całej macierzy
    %fft_matrix_current_step=fft_matrix_current_step/norm_factor;
    dlmwrite(current_fft_matrix_label, fft_matrix_current_step, 'newline', 'pc');
    
end
disp('Wykonano wszystkie transformaty i zapisano rezultaty do plików');
fprintf('Wykonano następującą liczbę transformat: %d\n', fft_done)

end


function draw_surface(pwr, Ts, Fs, ~, nmb_of_samples, g)
    for k=1:1:length(pwr)
        actual_fft_matrix_label=sprintf('fft_matrix/test/length_%d_plik_%d.txt', pwr(k),g);
        actual_fft_matrix = importdata(actual_fft_matrix_label, ',');
        
        fig_act_fft=surf(actual_fft_matrix,'EdgeColor','none');
        [row_nmb, col_nmb] = size(actual_fft_matrix);
        title_label=sprintf('Wykres macierzy transformaty\n o długości %d punktów - plik nr %d', col_nmb, g);
        title(title_label);
        fig_act_fft.XData=linspace(0,Fs,col_nmb);
        fig_act_fft.YData=linspace(0,Ts*nmb_of_samples,row_nmb);
        
        xlabel('f [Hz]');
        ylabel('t [s]');
        set(gca,'YDir','reverse');
        zlabel('widmo ampl.');
        fig_file_name=sprintf('fft_matrix/test/surface_%d_plik%d', pwr(k),g);
        saveas(fig_act_fft,fig_file_name, 'png');
        view([7, 70]); %view([0, 90]);
        fig_file_name=sprintf('fft_matrix/test/surface_%d_cv_plik%d', pwr(k),g);
        saveas(fig_act_fft,fig_file_name, 'png');
    end
end



