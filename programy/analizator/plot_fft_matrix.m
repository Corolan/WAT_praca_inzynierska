pwr = [8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072];   
Fs= 9765.625;
Ts=1/Fs;
nmb_of_samples=2^20;

k=10;
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
%fig_file_name=sprintf('fft_matrix/surface_%d', pwr(k));
%saveas(fig_act_fft,fig_file_name, 'fig');
%saveas(fig_act_fft,fig_file_name, 'png');
view([7, 70]);
%fig_file_name=sprintf('fft_matrix/surface_%d_cv', pwr(k));
%saveas(fig_act_fft,fig_file_name, 'png');