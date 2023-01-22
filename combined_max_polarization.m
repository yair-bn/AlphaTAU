function [snr] = combined_max_polarization(fig,window_size,NumOfSegment)
h = findobj(fig,'Type','line');
y = get(h,'Ydata');
x = get(h,'Xdata');
peaks = y{1};
peaks_loc = x{1};
noise = y{2};


% remove peaks:
for ii = [length(peaks):-1:1]

   cleaning_idx = max(floor(peaks_loc(ii)-window_size/2),1):min(ceil(peaks_loc(ii)+window_size/2), length(noise));  
   noise(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
 
end

% remove peak on the right side of the graph due to peak tail coming from the left edge: 
if peaks_loc(1)-window_size/2 < 1
    cleaning_idx = (length(noise)+(peaks_loc(1)-window_size/2)):length(noise);
    noise(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
end

% remove peak on the left side of the graph due to peak tail coming from the right edge: 
if peaks_loc(end)+window_size/2 > length(noise)
    cleaning_idx = 1: (peaks_loc(end)+window_size/2 - length(noise));
    noise(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
end

% REMOVE Inf values 
noise(noise==Inf)=[];


% convert to linear:
noise_lin = 10.^(noise./10);


% average noise (x,y):
avg_noise_lin = mean(noise_lin);


%convert to db:
avg_noise_dB = 10*log10(avg_noise_lin);


% calc square root of number of segments in dB:
SNR_improve = 10*log10(sqrt(NumOfSegment));

%signal - averaged noise:
snr = peaks - avg_noise_dB -SNR_improve;

figure; plot(noise,'linewidth',2); grid;
title(['noise level in dB', newline,num2str(avg_noise_dB),'dB'])
ylabel('noise [dB]'); xlabel('Samples')

end