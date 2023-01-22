function [xsnr,ysnr] = max_polarization(fig,window_size)
h = findobj(fig,'Type','line');
y = get(h,'Ydata');
x = get(h,'Xdata');
xpol_peaks = [x{3};y{3}];
ypol_peaks = [x{1};y{1}];

y_pol = [y{2};x{2}]; %[xpol(y);ypol(y);x]
x_pol = [y{4};x{4}]; %[xpol(y);ypol(y);x]

noise_x = x_pol(1,:);
noise_y = y_pol(1,:);

% remove peaks:
for ii = [length(x{3}):-1:1]

   cleaning_idx = max(xpol_peaks(1,ii)-window_size/2,1):min(xpol_peaks(1,ii)+window_size/2);length(noise_x);
   
   noise_x(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
   noise_y(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
   
   if xpol_peaks(1,ii)-window_size/2 < 1 
       cleaning_idx = length(noise_x)+(xpol_peaks(1,ii)-window_size/2):length(noise_x);
       
       noise_x(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
       noise_y(cleaning_idx) = ones(1,length(cleaning_idx))*Inf;
   end
end

% REMOVE Inf values 
noise_x(noise_x==Inf)=[];
noise_y(noise_y==Inf)=[];


% convert to linear:
noise_x_lin = 10.^(noise_x./10);
noise_y_lin = 10.^(noise_y./10);

% average noise (x,y):
avg_noise_lin = [mean(noise_x_lin),mean(noise_y_lin)];

%convert to db:
avg_noise_dB = 10*log10(avg_noise_lin);

%signal - averaged noise:
xsnr = xpol_peaks(2,:) - avg_noise_dB(1);
ysnr = ypol_peaks(2,:) - avg_noise_dB(2);


end