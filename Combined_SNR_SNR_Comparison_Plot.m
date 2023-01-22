clearvars
Folder = uigetdir('D:\alpha\Yair\exp_results\');
load([Folder,'\SNR'])
SNR_Leg_59 = combined_snr;
Folder = uigetdir('D:\alpha\Yair\exp_results\');
load([Folder,'\SNR'])
SNR_Leg_47 = combined_snr;
Folder = uigetdir('D:\alpha\Yair\exp_results\');
load([Folder,'\SNR'])
SNR_Ter_57 = combined_snr;

Folder = uigetdir('D:\alpha\Yair\exp_results\');

type = questdlg('Type of setup?','Setup', ...
    'Ring Fast','FBG Array', 'Ring Slow','Ring Fast');

if (strcmp(type,'Ring Fast'))
    [m,max_59] = max(SNR_Leg_59)
    [m,max_47] = max(SNR_Leg_47)
    [m,max_57] = max(SNR_Ter_57)
    
    peaks_order_Ter57 = [max_57:-1:1]; SNR_Ter_57=SNR_Ter_57(peaks_order_Ter57);
    peaks_order_Leg47 = [max_47:-1:1]; SNR_Leg_47=SNR_Leg_47(peaks_order_Leg47);
    peaks_order_Leg59 = [max_59:-1:1]; SNR_Leg_59=SNR_Leg_59(peaks_order_Leg59);
    smallest_code = min(length(peaks_order_Ter57), length(peaks_order_Leg47))
    smallest_code = min(smallest_code, length(peaks_order_Leg59))
end

if (strcmp(type,'Ring Slow'))
    smallest_code = min(length(SNR_Ter_57), length(SNR_Leg_47))
    smallest_code = min(smallest_code, length(SNR_Leg_59))
end

MeanLeg59_AllSNR = mean(Leg59_lin(:));
MeanLeg47_AllSNR = mean(Leg47_lin(:));
MeanTer57_AllSNR = mean(Ter57_lin(:));
    
disp('Ter57 vs. Leg59 (All values)')
disp(MeanTer57_AllSNR - MeanLeg59_AllSNR)
disp('Ter57 vs. Leg47 (All values)')
disp(MeanTer57_AllSNR - MeanLeg47_AllSNR)
    
figure(1);
plot(SNR_Ter_57,'c*','linewidth',2.5); grid
hold on; plot(SNR_Leg_59,'g*','linewidth',2.5);
hold on; plot(SNR_Leg_47,'k*','linewidth',2.5);
title('SNR Comparison for Ternary 57, Legendre 59 and Legendre 47','fontsize',20)
xlabel('Reflector number','fontsize',18)
ylabel('Static SNR (Per Segment)','fontsize',18)
if (strcmp(type,'Ring Fast') || strcmp(type,'Ring Slow'))
    xlim([0 (smallest_code+0.9)])
end
legend('Ternary 57','Legendre 59','Legendre 47')
savefig(1,[Folder,'\SNR_Comparison_Ternary_57_Legendre_59_Legendre_47'])

close all;