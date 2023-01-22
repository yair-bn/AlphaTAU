clearvars
Folder = uigetdir('D:\alpha\Yair\FINAL_RESULTS\Experiment');
load([Folder,'\SNR'])
SNR_Leg_59 = combined_snr;
Folder = uigetdir('D:\alpha\Yair\FINAL_RESULTS\Experiment');
load([Folder,'\SNR'])
SNR_Leg_47 = combined_snr;
Folder = uigetdir('D:\alpha\Yair\FINAL_RESULTS\Experiment');
load([Folder,'\SNR'])
SNR_Ter_57 = combined_snr;

sum_SNR_Leg_59 = 0;
sum_SNR_Leg_47 = 0;
sum_SNR_Ter_57 = 0;



for i = 1:length(SNR_Leg_59)
    sum_SNR_Leg_59 = sum_SNR_Leg_59 + SNR_Leg_59(i);
end

for i = 1:length(SNR_Leg_47)
    sum_SNR_Leg_47 = sum_SNR_Leg_47 + SNR_Leg_47(i);
end

for i = 1:length(SNR_Ter_57)
    sum_SNR_Ter_57 = sum_SNR_Ter_57 + SNR_Ter_57(i);
end


avg_SNR_Leg_59 = sum_SNR_Leg_59 / length(SNR_Leg_59);
avg_SNR_Leg_47 = sum_SNR_Leg_47 / length(SNR_Leg_47);
avg_SNR_Ter_57 = sum_SNR_Ter_57 / length(SNR_Ter_57);

dif_SNR_Leg_59_Ter_57 = avg_SNR_Ter_57 - avg_SNR_Leg_59;
dif_SNR_Leg_47_Ter_57 = avg_SNR_Ter_57 - avg_SNR_Leg_47;

disp(dif_SNR_Leg_59_Ter_57)
disp(dif_SNR_Leg_47_Ter_57)