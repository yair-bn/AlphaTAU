clearvars
close all;

prompt = {'Enter number of Leg 59 experiments:'};
definput = {'1'}; dims = [1 44]; opts.Resize = 'on';
aa = inputdlg2(prompt,'Experiments Num',dims,definput,opts);
NumOfExpLeg59 = str2num(aa{1}); 
clear aa;

leg59_idx = {};
for ii=1:NumOfExpLeg59
    if ii==1
        Folder = uigetdir('D:\alpha\Yair\exp_results\');
        temp = find(Folder=='\');
        Folder_focus = Folder(1:temp(end));
    else
        Folder = uigetdir(Folder_focus);
    end
    load([Folder,'\SNR'])
    SNR_Leg_59 = combined_snr;
    if ii==1
        smallestLeg59 = length(SNR_Leg_59);
    end
    if (length(SNR_Leg_59)<smallestLeg59)
        smallestLeg59=length(SNR_Leg_59);
    end
    leg59_idx{end+1} = SNR_Leg_59;
end



prompt = {'Enter number of Leg 47 experiments:'};
definput = {'1'}; dims = [1 44]; opts.Resize = 'on';
aa = inputdlg2(prompt,'Experiments Num',dims,definput,opts);
NumOfExpLeg47 = str2num(aa{1}); % in meters
clear aa;

leg47_idx = {};

for ii=1:NumOfExpLeg47
  
    Folder = uigetdir(Folder_focus);    
    load([Folder,'\SNR'])
    SNR_Leg_47 = combined_snr;
    if ii==1
        smallestLeg47 = length(SNR_Leg_47);
    end
    if (length(SNR_Leg_47)<smallestLeg47)
        smallestLeg47=length(SNR_Leg_47);
    end
    leg47_idx{end+1} = SNR_Leg_47;
end



prompt = {'Enter number of Ter 57 experiments:'};
definput = {'1'}; dims = [1 44]; opts.Resize = 'on';
aa = inputdlg2(prompt,'Experiments Num',dims,definput,opts);
NumOfExpTer57 = str2num(aa{1}); % in meters
clear aa;


ter57_idx = {};
for ii=1:NumOfExpTer57
    
    Folder = uigetdir(Folder_focus);
    load([Folder,'\SNR'])
    SNR_Ter_57 = combined_snr;
    if ii==1
        smallestTer57 = length(SNR_Ter_57);
    end
    if (length(SNR_Ter_57)<smallestTer57)
        smallestTer57=length(SNR_Ter_57);
    end
    ter57_idx{end+1} = SNR_Ter_57;
end

Folder = uigetdir('D:\alpha\Yair\exp_results\', 'Enter Folder For Saving');

type = questdlg('Type of setup?','Setup', ...
    'Ring Fast' , 'FBG Array', 'Ring Slow','Ring Fast');

if (strcmp(type,'Ring Fast'))
    peaks_order_Leg59_idx = {}
    peaks_order_Leg47_idx = {}
    peaks_order_Ter57_idx = {}
    
%     for ii=1:NumOfExpLeg59
%         [m,max_59] = max(leg59_idx{ii})
%         peaks_order_Leg59_idx{end+1} = {max_59:-1:1};leg59_idx{ii}=leg59_idx{ii(peaks_order_Leg59_idx{ii})};
%     end 


    [m,max_59] = max(SNR_Leg_59)
    peaks_order_Leg59 = [max_59:-1:1]; SNR_Leg_59=SNR_Leg_59(peaks_order_Leg59);
    
    [m,max_47] = max(SNR_Leg_47);
    peaks_order_Leg47 = [max_47:-1:1]; SNR_Leg_47=SNR_Leg_47(peaks_order_Leg47);
    
    [m,max_57] = max(SNR_Ter_57);
    peaks_order_Ter57 = [max_57:-1:1]; SNR_Ter_57=SNR_Ter_57(peaks_order_Ter57);
    
    smallest_code = min(length(peaks_order_Ter57), length(peaks_order_Leg47))
    smallest_code = min(smallest_code, length(peaks_order_Leg59))
end

if (strcmp(type,'Ring Slow'))
    smallest_code = min( min(smallestTer57, smallestLeg59), min(smallestTer57, smallestLeg47) );
end

if (strcmp(type,'FBG Array'))
    smallest_code = length(SNR_Ter_57);
end
%% calc average SNR for each code

% from cell array to matrix:
for kk = 1:NumOfExpLeg59
    Leg59_dB(kk,:) = leg59_idx{kk}(1:smallest_code);
end

for kk = 1:NumOfExpLeg47
    Leg47_dB(kk,:) = leg47_idx{kk}(1:smallest_code);
end

for kk = 1:NumOfExpTer57
    Ter57_dB(kk,:) = ter57_idx{kk}(1:smallest_code);
end

% from dB to linear:
Leg59_lin = 10.^(Leg59_dB/10);
Leg47_lin = 10.^(Leg47_dB/10);
Ter57_lin = 10.^(Ter57_dB/10);

% calc average and std :
MeanLeg59 = mean(Leg59_lin);  StdLeg59 = std(Leg59_lin);
MeanLeg47 = mean(Leg47_lin);  StdLeg47 = std(Leg47_lin);
MeanTer57 = mean(Ter57_lin);  StdTer57 = std(Ter57_lin);


% calc average SNR for SNR > 3 :
if (strcmp(type,'Ring Fast')||strcmp(type,'Ring Slow'))
    %big_3 = gt(MeanLeg47,2)
    %NumOfRef = sum(big_3(:) == 1)
    NumOfRef = 14
    MeanLeg59_HighSNR = mean(MeanLeg59(1:NumOfRef));  %StdLeg59_HighSNR = std(Leg59_lin(1:NumOfRef));
    MeanLeg47_HighSNR = mean(MeanLeg47(1:NumOfRef));  %StdLeg47_HighSNR = std(Leg47_lin(1:NumOfRef));
    MeanTer57_HighSNR = mean(MeanTer57(1:NumOfRef));  %StdTer57_HighSNR = std(Ter57_lin(1:NumOfRef));
    
    MeanLeg59_LowSNR = mean(MeanLeg59(NumOfRef+1:end));
    MeanLeg47_LowSNR = mean(MeanLeg47(NumOfRef+1:end));
    MeanTer57_LowSNR = mean(MeanTer57(NumOfRef+1:end)); 
    
    MeanLeg59_AllSNR = mean(MeanLeg59(:));
    MeanLeg47_AllSNR = mean(MeanLeg47(:));
    MeanTer57_AllSNR = mean(MeanTer57(:));
    
    disp('Ter57 vs. Leg59 (High values)')
    disp(10*log10(MeanTer57_HighSNR / MeanLeg59_HighSNR))
    disp('Ter57 vs. Leg47 (High values)')
    disp(10*log10(MeanTer57_HighSNR / MeanLeg47_HighSNR))
    
    disp('Ter57 vs. Leg59 (Low values)')
    disp(10*log10(MeanTer57_LowSNR / MeanLeg59_LowSNR))
    disp('Ter57 vs. Leg47 (Low values)')
    disp(10*log10(MeanTer57_LowSNR / MeanLeg47_LowSNR))
    
    disp('Ter57 vs. Leg59 (All values)')
    disp(10*log10(MeanTer57_AllSNR / MeanLeg59_AllSNR))
    disp('Ter57 vs. Leg47 (All values)')
    disp(10*log10(MeanTer57_AllSNR / MeanLeg47_AllSNR))
else
    MeanLeg59_AllSNR = mean(MeanLeg59(:));
    MeanLeg47_AllSNR = mean(MeanLeg47(:));
    MeanTer57_AllSNR = mean(MeanTer57(:));
    
    disp('Ter57 vs. Leg59 (All values)')
    disp(10*log10(MeanTer57_AllSNR / MeanLeg59_AllSNR))
    disp('Ter57 vs. Leg47 (All values)')
    disp(10*log10(MeanTer57_AllSNR / MeanLeg47_AllSNR))
end
%% my plots

type = questdlg('Type of comparison?','comparison', ...
    'Same transmitted energy' , 'Same T Scan','Same T Scan');
% linear plot:
figure(12)
h1 = errorbar(MeanLeg59,StdLeg59,'o'); hold on;
h2 = errorbar(MeanLeg47,StdLeg47,'o'); hold on;
h3 = errorbar(MeanTer57,StdTer57,'o'); grid on;
legend('Legendre 59','Legendre 47','Ternary 57')
h1.LineWidth=2; h2.LineWidth=2;  h3.LineWidth=2; 
ylabel('SNR [Linear]'); xlabel('Number of rounds');
title(['SNR Comparison - ',type],'fontsize',20)
ax=gca; ax.LineWidth=1.5;
ax.FontSize=16;

% log plot:
figure(13)
h1 = plot(10*log10(MeanLeg59),'o'); hold on;
h2 = plot(10*log10(MeanLeg47),'o'); hold on;
h3 = plot(10*log10(MeanTer57),'o'); grid on;
legend('Legendre 59','Legendre 47','Ternary 57');
h1.LineWidth=2; h2.LineWidth=2;  h3.LineWidth=2; 
ylabel('SNR [dB]'); xlabel('Number of rounds');
title(['SNR Comparison - ',type],'fontsize',20)
ax=gca; ax.LineWidth=1.5; ax.FontSize=16;

% save figure:
savefig(12,[Folder,'\SNR_Comparison_Ternary_57_Legendre_59_Legendre_47_Linear',type(end-6:end)])
savefig(13,[Folder,'\SNR_Comparison_Ternary_57_Legendre_59_Legendre_47_dB',type(end-6:end)])
%% Yair's plot
figure(1);
Legend = {};
for ii=1:NumOfExpTer57
    plot(ter57_idx{ii},'c*','linewidth',2.5); grid
    hold on  
    Legend{end+1} = 'Ternary 57';
end
for ii=1:NumOfExpLeg59
    plot(leg59_idx{ii},'g*','linewidth',2.5);
    hold on
    Legend{end+1} = 'Legendre 59';
end
for ii=1:NumOfExpLeg47
    plot(leg47_idx{ii},'k*','linewidth',2.5);
    hold on; 
    Legend{end+1} = 'Legendre 47';
end
title('SNR Comparison for Ternary 57, Legendre 59 and Legendre 47','fontsize',20)
xlabel('Reflector number','fontsize',18); grid on;
ylabel('Static SNR (Per Segment)','fontsize',18)
if (strcmp(type,'Ring Fast') || strcmp(type,'Ring Slow'))
    xlim([0 (smallest_code+0.9)])
end
legend(Legend)

