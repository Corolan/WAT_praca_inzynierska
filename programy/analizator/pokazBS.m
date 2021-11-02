tic;
clear variables; 
filename = sprintf('wynik.txt');
raw_data = importdata(filename,' ');
Re = raw_data(:,1);
Im = raw_data(:,2);
data=complex(Re, Im);

data_length=length(data);
bsmooth_boxcar_length=8192;
bsmooth_chunk_size=32768*32;
Fs= 9765.625;

%f=linspace(0,Fs,length(data));
length(data)
f=linspace(0,2^20,length(data));
xl='częstotliwość [Hz]';
yl='Moduł transformaty fouriera';
    
num_of_smooths = data_length/bsmooth_chunk_size;
frequency_domain_before_smoothing = (fft(data,data_length)).^2;
figure;
subplot(1,2,1);
przedBS=plot(f,abs(frequency_domain_before_smoothing));
%przedBS=plot(abs(data), 'r*');
axis([-inf inf -inf inf]);
title('Dane przed wygładzeniem');
ylabel(yl);
xlabel(xl);
  
for cnt=1:1:num_of_smooths
    act_data=frequency_domain_before_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt); %aktualnie przetwarzany kawałek danych
    frequency_domain_after_smoothing((bsmooth_chunk_size)*(cnt-1)+1:(bsmooth_chunk_size)*cnt)=(smooth(act_data, bsmooth_boxcar_length))';
end
    
subplot(1,2,2);
data_after=ifft(frequency_domain_after_smoothing, length(frequency_domain_after_smoothing));
%poBS=plot(f,abs(frequency_domain_after_smoothing));
data_to_plot=data_after;

%data_to_plot=data_after(200000:length(data_after)-200000);

poBS=plot(abs(data_to_plot), 'g*');
axis([-inf inf -inf inf]);
title('Dane po wygładzeniu');
ylabel(yl);
xlabel(xl);

% subplot(1,3,3);
% poBS1=plot(f,abs(frequency_domain_after_smoothing)/max(abs(frequency_domain_after_smoothing)));
% %poBS=plot(abs(frequency_domain_after_smoothing)/max(abs(frequency_domain_after_smoothing)));
% axis([-inf inf 0 0.1]);
% title_label=sprintf('Dane po wygładzeniu\n (zmiana skali osi Y)');
% title(title_label);
% ylabel(yl);
% xlabel(xl);

saveas(gcf, 'bs_smooth/baseline_smoothing-show', 'jpeg');
toc;
