% main <not in real time >
clc
close all
clearvars
debug = 0;

%% enter folder location
go_back = 1;
while(go_back == 1) % keep looping if there is no h5 file in the folder entered by the user
    
    Folder = uigetdir('D:\Nadav_D\experiments\FBG_array');
    file = dir([Folder ,'/*h5']);
    filename = [file(1).folder ,'\',file(1).name];
    
    if ~isempty(file)
        go_back=0;
        if size(file,1)>1
            msgbox('there is more then one h5 file in the folder')
        end
    else
        disp('there is no h5 file format inside the folder')
    end
    
end

%% chose reference for the cross-correlation process

aa = questdlg('\fontsize{15}type of reference?','XCORR', ...
    'Digital','External',struct('Default','Digital','Interpreter','tex'));
if strcmp(aa,'Digital')
    use_external_waveform = 0;
else
    use_external_waveform = 1;
end
clear aa;

complex_ref =0;

%% set params

SNRTH = 10;        % min static SNR to use

%% corr factor -> for IQ imbalance

corr_factor_x = 1; %681/927;
corr_factor_y = 1; %1380/1190;

X_corr_flag = 1;
Y_corr_flag = 1;

if corr_factor_x <1
    corr_factor_x = 1/corr_factor_x;
    X_corr_flag = 1;
end
if corr_factor_y <1
    corr_factor_y = 1/corr_factor_y;
    Y_corr_flag = 1;
end

%% create a result subfolder
if use_external_waveform==1
    results_folder = [file(1).folder,'\results_ext_ref'];
elseif use_external_waveform==0
    results_folder = [file(1).folder,'\results_dig_ref'];
end

if ~exist(results_folder, 'dir')
    mkdir(results_folder)
end
%% scripts parameters

NumOfChan = double(h5readatt(filename,'/Waveforms','NumWaveforms'));
I=1;
Q=I;

if NumOfChan<4
    
    while I==Q
        gen_data = h5info(filename);
        gen_data.Groups(3).Groups.Name
        
        for nn=1:NumOfChan
            temp = gen_data.Groups(3).Groups(nn).Name;
            opt(nn) = temp(end);
        end
        
        if NumOfChan==3
            I = questdlg('Enter I channel ', ...
                'Enter I channel', ...
                opt(1),opt(2),opt(3),opt(1));
            
            
            Q = questdlg('Enter Q channel ', ...
                ' I Q', ...
                opt(1),opt(2),opt(3),opt(1));
        else
            I = questdlg('Enter I channel ', ...
                'Enter I channel', ...
                opt(1),opt(2),opt(1));
            
            
            Q = questdlg('Enter Q channel ', ...
                ' I Q', ...
                opt(1),opt(2),opt(1));
        end
    end
    
elseif NumOfChan == 4
    
    I=1;
    Q=2;
    
end

aa = questdlg('\fontsize{15}Enter array roundtrip time ', ...
    'FBG', ...
    '490.2','195.8',struct('Default','195.8','Interpreter','tex'));
T_RoundTrip = str2num(aa);
clear aa

aa = questdlg('\fontsize{15}Enter Legendre length ', ...
    'Legendre', ...
    '83','199','other',struct('Default','199','Interpreter','tex'));
code_length = str2num(aa);
clear aa


if isempty(code_length)
    
    code_length = 0;
    while ~isprime(code_length)
        prompt = {'\fontsize{15}Enter code length [must be a prime number]'};
        definput = {'103'}; dims = [1 44]; opts.Resize = 'on';
        aa = inputdlg2(prompt,'Legendre bits',dims,definput,opts);
        code_length = str2num(aa{1}); % in meters
        clear aa
    end
    
end

est_pulse_dur = T_RoundTrip/code_length; % [ns]
est_pulse_dur_round = round(est_pulse_dur*10)/10; % [ns]

if T_RoundTrip == 490.2
    
    % for 50m
    if (est_pulse_dur >= 2) && (est_pulse_dur <= 3)
        opt = {'2.5','2','other'};
    elseif (est_pulse_dur >= 5) && (est_pulse_dur <= 6)
        opt = {'5.9','6','other'};
    else
        opt = {num2str(est_pulse_dur_round),num2str(est_pulse_dur_round+0.1),'other'};
    end
    
else
    % for 20m
    
    if code_length==199
        pulse_dur = 1;   % [ns]
    elseif code_length==83
        pulse_dur = 2.4; % [ns]
    else
        opt = {num2str(est_pulse_dur_round),num2str(est_pulse_dur_round+0.1),'other'};
    end
    
end

if ~exist('pulse_dur')
    aa = questdlg('\fontsize{15}Enter pulse duration in ns', ...
        'Legendre', ...
        opt{1},opt{2},opt{3},struct('Default',opt{1},'Interpreter','tex'));
    pulse_dur = str2num(aa);
    clear aa opt
end


if isempty(pulse_dur)
    prompt = {'\fontsize{15}Enter pulse duration in ns'};
    definput = {num2str(est_pulse_dur_round)}; dims = [1 40]; opts.Resize = 'on';
    aa = inputdlg2(prompt,'Legendre bits',dims,definput,opts);
    pulse_dur = str2num(aa{1}); % in nano-sec
    clear aa
end


if NumOfChan <4
    
    if (str2num(I)==1)&&(str2num(Q)==2)
        aa = ['3','4','4'];
    elseif (str2num(I)==1)&&(str2num(Q)==3)
        aa = ['2','4','4'];
    elseif (str2num(I)==1)&&(str2num(Q)==4)
        aa = ['2','3','3'];
    elseif (str2num(I)==2)&&(str2num(Q)==3)
        aa = ['1','4','4'];
    end
    
    bb = questdlg('\fontsize{15}Enter trigger channel number ', ...
        'Legendre', ...
        aa(1),aa(2),'not exist',struct('Default',aa(3),'Interpreter','tex'));

    T = str2num(bb);
    clear aa bb
    
end

pulse_dur;                       % [nsec]
pulse_duration = pulse_dur*1e-9; % [sec]

%% load recorded data of .h5 file type
disp('phase2 - load raw data')

I_data_x = double(h5read( filename, ['/Waveforms/Channel ',num2str(I),'/Channel ',num2str(I),'Data' ]));
Q_data_x = double(h5read( filename, ['/Waveforms/Channel ',num2str(Q),'/Channel ',num2str(Q),'Data' ]));

if NumOfChan==4
    I_data_y = double(h5read( filename, ['/Waveforms/Channel ',num2str(3),'/Channel ',num2str(3),'Data' ]));
    Q_data_y = double(h5read( filename, ['/Waveforms/Channel ',num2str(4),'/Channel ',num2str(4),'Data' ]));
    I_data_y = I_data_y - mean(I_data_y); % subtract DC component from raw data
    Q_data_y = Q_data_y - mean(Q_data_y); % subtract DC component from raw data
end

% subtract DC component from raw data
% only for coding interrogation as there is no DC
I_data_x = I_data_x - mean(I_data_x);
Q_data_x = Q_data_x - mean(Q_data_x);

%% fast time axis
Fs = 1/double(h5readatt(filename,['/Waveforms/Channel ',num2str(I)],'XInc'));  % sampling rate
dt = 1/Fs;                           % Sampling time duration [sec]
t = (0:dt:(length(I_data_x)-1)*dt);  % [sec]

SegLength = double(h5readatt(filename,['/Waveforms/Channel ',num2str(I)],'NumPoints'));% num of samples

%%
if (debug)
    
    figure; plot(t*1e9,I_data_x,'b','linewidth',2.5); grid
    hold on; plot(t*1e9,Q_data_x,'r','linewidth',2.5);
    xlabel('time [ns]','fontsize',16)
    xlim([t(1) t(end)]*1e9)
    legend('I X pol','Q X pol','I Y pol','Q Y pol')
    
    figure;  plot(t(1:5e3)*1e6,I_data_x(1:5e3),'b','linewidth',2.5); grid
    hold on; plot(t(1:5e3)*1e6,Q_data_x(1:5e3),'r','linewidth',2.5);
    xlabel('time [\mus]','fontsize',16)
    xlim([t(1) t(5e3)]*1e6)
    legend('I X pol','Q X pol')
    
    data_amp = abs(I_data_x + 1i*Q_data_x);
    data_p = atan(I_data_x./Q_data_x);
    
    figure;
    plot(t*1e6,data_amp,'b','linewidth',2.5); grid
    title('raw data','fontsize',16)
    xlabel('time [us]','fontsize',16)
    xlim([t(1) t(end)]*1e6)
    
    figure;
    plot(t*1e6,data_p,'b','linewidth',2.5); grid
    title('phase of raw data','fontsize',16)
    xlabel('time [us]','fontsize',16)
    ylabel('phase [rad]','fontsize',16)
    xlim([t(1) t(end)]*1e6)
end

%% correct for imbalnce between I and Q
disp('phase3 - correct for IQ unbalance')

% [corr_factor]  = IQbalance(I_data_x,Q_data_x);  % [corr_factor2] = IQbalance(I_data_x,Q_data_x);

clear  I_data_x_temp Q_data_x_temp
TH = max(abs(I_data_x+1i*Q_data_x));
I_data_x_temp = I_data_x (abs(I_data_x+1i*Q_data_x)>0.65*TH);
Q_data_x_temp = Q_data_x (abs(I_data_x+1i*Q_data_x)>0.65*TH);

if X_corr_flag==0
    figure; scatter (I_data_x_temp ,corr_factor_x*Q_data_x_temp )
else
    figure; scatter (corr_factor_x*I_data_x_temp ,Q_data_x_temp )
end

grid; text(0,0,['corr factor = ' ,num2str(corr_factor_x)]);
xlim([min(I_data_x_temp) , max(I_data_x_temp)]);
ylim([min(I_data_x_temp) , max(I_data_x_temp)]);% axis('equal');
if X_corr_flag == 0
    data_c_x = I_data_x + corr_factor_x*1i*Q_data_x;
else
    data_c_x = corr_factor_x*I_data_x + 1i*Q_data_x;
end

if NumOfChan ==4
    
    TH = max(abs(I_data_y+1i*Q_data_y));
    I_data_y_temp = I_data_y (abs(I_data_y+1i*Q_data_y)>0.8*TH);
    Q_data_y_temp = Q_data_y (abs(I_data_y+1i*Q_data_y)>0.8*TH);
    
    if Y_corr_flag==0
        figure; scatter (I_data_y_temp ,corr_factor_y*Q_data_y_temp )
    else
        figure; scatter (corr_factor_y*I_data_y_temp ,Q_data_y_temp )
    end
    grid; text(0,0,['corr factor = ' ,num2str(corr_factor_y)])
    xlim([min(I_data_y_temp) , max(I_data_y_temp)]);
    ylim([min(I_data_y_temp) , max(I_data_y_temp)]);% axis('equal');
    if Y_corr_flag==0
        data_c_y = I_data_y + corr_factor_y*1i*Q_data_y;
    else
        data_c_y = corr_factor_y*I_data_y + 1i*Q_data_y;
    end
end

clear  I_data_x_temp Q_data_x_temp
close all;

%% cross-correlation of raw data with reference Legendre seq
disp('phase4 - cross correlation')

if use_external_waveform==0
    
    [waveform] = gen_Legendre_seq(code_length,pulse_dur,Fs);
    NumOfCopies = 1; % note it effects acoustic BW
    waveform_ex = repmat(waveform, 1, NumOfCopies);
    
elseif use_external_waveform==1
    go_back = 1;
    while(go_back == 1) % keep looping if there is no h5 file in the folder entered by the user
        
        if T_RoundTrip==490.2
            Folder = uigetdir('D:\Alon_David\experiments_results\Experiment_Results\refrence_signal\FBG50');
        else
            Folder = uigetdir('D:\Alon_David\experiments_results\Experiment_Results\refrence_signal\FBG20');
        end
        file = dir([Folder ,'/*h5']);
        reffilename = [file(1).folder ,'\',file(1).name];
        
        if ~isempty(file)
            go_back=0;
            if size(file,1)>1
                msgbox('there is more then 1 h5 file in the folder')
            end
        else
            disp('there is no h5 file format inside the folder')
        end
    end
    waveform_ex_I = double(h5read( reffilename, ['/Waveforms/Channel ',num2str(1),'/Channel ',num2str(1),'Data' ]));
    if complex_ref == 1
        waveform_ex_Q = double(h5read( reffilename, ['/Waveforms/Channel ',num2str(2),'/Channel ',num2str(2),'Data' ]));
        waveform_ex = waveform_ex_I + 1i*waveform_ex_Q;
    else
        waveform_ex = waveform_ex_I;
    end
    NumOfCopies = 1; % note it effect acoustic BW
    
end


% option A :
detect_x = xcorr(data_c_x ,waveform_ex );   %
detect_x = detect_x(round(end/2-0.5*length(waveform_ex)) : end - round(0.5*length(waveform_ex)));

% option B :
% detect_x2 = cconv(data_c_x,conj(fliplr(waveform_ex)),length(data_c_x));

% option C :
% detect_x = ifft(fft(data_c_x) .* conj(fft(waveform)));

detect_db = dB20(detect_x);

if NumOfChan ==4
    detect_y = xcorr(data_c_y ,waveform_ex );
    detect_y = detect_y(round(end/2-0.5*length(waveform_ex)) : end - round(0.5*length(waveform_ex)));
end

%% find zero-crossing in order to fold the vector into a waterfall image

period = round(code_length*pulse_duration*Fs); % [samples]
zcros = zeros(1,round(length(detect_x)/period)-2); space = 45;

for ii=1:length(zcros)-1
    
    if ii==1
        
        
        findpeaks(abs(detect_x( period + 1 : 2*period+5 )),'MinPeakHeight' ,max(abs(detect_x(period+1:2*period)))/5);
        title('Set Y Threshold for peak detection');temp = ginput(1); close;  Y_TH_lin = temp(2); clear temp;
        [peak_height,pos] = findpeaks(abs(detect_x( period + 1 : 2*period+5 )),'MinPeakHeight' , Y_TH_lin);
        pos(pos>period)=[]; % clear peaks that are outside the period
        
        % find the highest peak
        [a1,peak_idx] = max(peak_height);
        pos1 = pos(peak_idx);
        
        
        if (length(pos))== 1
            diff_vec = [pos(1)]; % H-spacing between peaks [samps]
        else
            diff_vec = diff(pos);
        end
        space = round(min(diff_vec)/2);
        
        
        zz(ii) = pos1;
        zcros(ii) = pos1 + period;
        
        
    else
        
        
        % [~,pos] = max(abs(detect_x( 1+ period*(i) + pos1 - space : period*(i) + pos1  + space )));
        [val_peak,pos_list] = findpeaks(abs(detect_x( 1+ period*(ii) + pos1 - space : period*(ii) + pos1  + space ))) ;
        [~,max_peak_idx] = max(val_peak);
        pos = pos_list(max_peak_idx);
        
        zz(ii) = pos1 + pos - space;
        zcros(ii) = pos1 + pos - space + period*(ii);
        
        % update pos1 value:
        pos1 = mod(zz(ii),period);
    end
    
end

%%% find anomalous in peaks position %%%

% peak that moved one step to the right will have zz value higher than its neighbors:
[zz_val,zz_anom_idx] = findpeaks(zz,'MaxPeakWidth',2);
zz(zz_anom_idx) = zz(zz_anom_idx)-1;
zcros(zz_anom_idx) = zcros(zz_anom_idx)-1;

% peak that moved one step to the right will have zz value higher than its neighbors:
[zz_val_neg,zz_anom_idx_neg] = findpeaks(-zz,'MaxPeakWidth',2);
zz(zz_anom_idx_neg) = zz(zz_anom_idx_neg)+1;
zcros(zz_anom_idx_neg) = zcros(zz_anom_idx_neg)+1;

while(zcros(end)==0)
    zcros(end)=[];
end

if debug
    figure; plot(zz,'*');
    hold on; plot(zz_anom_idx,zz_val,'ro')
    hold on; plot(zz_anom_idx_neg,-zz_val_neg,'ko')
    title('peak position rel to segment starting point')
    
    figure
    plot(abs(detect_x(1:round(end/30)))); grid; hold on;
    plot(zcros(1:round(end/30)),abs(detect_x(zcros(1:round(end/30)))),'*')
end

% move the zero-cross point away from the peak, so the peak will be centered

figure;
plot(abs(detect_x(1:round(end/40)))); grid; hold on;
plot(zcros(1:round(end/40)),abs(detect_x(zcros(1:round(end/40)))),'*')
xlim([ round(length(detect_x)/(2*40))-round(period) , round(length(detect_x)/(2*40))+round(period)])
title('Choose a point between the two segments')
user_corr = ginput(1); horiz_pos = user_corr(1); close;

temp = find (zcros<horiz_pos,1,'last');
zcros_chosen = zcros(temp);
zcros_chosen(zcros_chosen==0)=[];

horiz_shift = round(horiz_pos - zcros_chosen);
zcros = zcros + horiz_shift;

%% calc min distance between two adjacent reflections
% code & bit lengths
bit_L = round(pulse_duration*Fs);   % in samples
period = bit_L*code_length;         % in sampels
first_peaks = 3;
samples = [];

min_dis_samp = min( (abs(code_length*pulse_dur-T_RoundTrip)*1e-9 )*Fs , median(diff_vec));  % in samples
min_dis_samp = min_dis_samp-11;                              % for tolerance only

%% reshape data to waterfall
% period =  median(diff(zcros));  % in samples

c_data_x = zeros ( length(zcros)-1, period );% Num Of codes recorded X code length in samples

if exist('detect_y','var')
    
    c_data_y = zeros ( length(zcros)-1, period );% Num Of codes recorded X code length in samples
    
    for kk = 1:(length(zcros)-1)
        
        c_data_x(kk,:) = detect_x(zcros(kk): zcros(kk) + period -1);
        c_data_y(kk,:) = detect_y(zcros(kk): zcros(kk) + period -1);
        
    end
    
else
    
    for kk = 1:(length(zcros)-1)
        
        c_data_x(kk,:) = detect_x(zcros(kk): zcros(kk) + period -1);
        
    end
    
end
c_data_db_x = dB20(c_data_x);

%% detected data paramaters
disp('phase5 - detected data paramaters')

% fast axis params :

% Fs;                 % Sampling frequency
% dt = 1/Fs;          % Sampling time duration (sec)
L = size(c_data_x,2); % number of samples per scan
T = dt*L;             % total time duration of a single scan / single sequence (sec)
t_fast = 0:dt:T-dt;
% freq parameters in the fast axis:
% df_fast = Fs/L;            % freq resolution (Hz)
% f_fast = -Fs/2 : df_fast : (Fs/2-df_fast); % freq vector (Hz) =  Fs*(-(L/2):1:(L/2)-1)/L;

% slow axis params :

T; % Sampling time duration [sec]
NumOfCodePeriods = floor(length(detect_db)/period);% [number of scans]-> in Theory
scan_num = size(c_data_x,1);   % [number of scans]-> in Practice
T_tot = T*scan_num;            % [sec]
t_slow = 0:T:T_tot-T;
df_slow = 1/T_tot;
Fs_slow = 1/T;
dt_slow = 1/Fs_slow;
f_slow = -Fs_slow/2 : df_slow : (Fs_slow/2-df_slow); % freq vector (Hz)

%% plot waterfall in time domain
if (debug)
    
    % plot waterfall in time domain
    figure(99)
    imagesc (t_fast,t_slow*1e3,c_data_db_x - max(c_data_db_x(:)),[ mean(c_data_db_x(:))-max(c_data_db_x(:))-10, mean(c_data_db_x(:))-max(c_data_db_x(:))+10] );
    colorbar
    ylabel('slow time [ms]','fontsize',16)
    xlabel('fast time (mod)','fontsize',16)
    title('abs of detected data in dB ','fontsize',18)
    
    % single segment in time domain
    figure(9); plot(t_fast,abs(detect_x(zcros(kk): zcros(kk) + median(diff(zcros))-1))); grid on;
    title('single segment','fontsize',18);ylabel('linear units','fontsize',16)
    xlabel('fast time in sec (mod)','fontsize',16); xlim([0,t_fast(end)])
    
    % single segment in time domain
    figure(19); plot(t_fast*2e8,abs(detect_x(zcros(kk): zcros(kk) + median(diff(zcros))-1))); grid on;
    title('single segment','fontsize',18); ylabel('linear units','fontsize',16)
    xlabel('distance in meters (mod)','fontsize',16); xlim([0,t_fast(end)]*2e8);
    
    % single segment in time domain -> linear and log scale together
    figure(29); yyaxis left;
    plot(t_fast*2e8,abs(c_data_x(end,:)),'-','linewidth',2.5); grid on;
    xlabel('distance in meters (mod)','fontsize',16); xlim([0,t_fast(end)]*2e8);
    ylabel('linear units','fontsize',16)
    yyaxis right;
    plot(t_fast*2e8,abs(c_data_db_x(end,:)),'--');
    title('single segment','fontsize',18);
    ylabel('dB','fontsize',16)
    
end

%% fft on complex data
if (debug)
    
    data_f_dsb = (fftshift(fft(c_data_x),1));        % two-sided spectrum
    data_f_dsb = data_f_dsb / max(abs(data_f_dsb(:)));  % normalization
    
    
    figure(100); imagesc( t_fast,f_slow,20*log10(abs(data_f_dsb)) )
    colorbar
    ylabel('Freq [Hz]','fontsize',16)
    xlabel('Fast time (mod)','fontsize',16)
    title('Fourier transform of detected data [dB]','fontsize',18)
end

%% calc the average power of all segments
disp('phase6 - calc average power of all segments')
% average over power:
mean_seg_x = mean(abs(c_data_x).^2,1);      % mean_seg_x = mean(abs(c_data_x(1:250,:)).^2,1);

if NumOfChan == 4
    % average over power
    mean_seg_y = mean(abs(c_data_y).^2,1);  % mean_seg_y = mean(abs(c_data_y(1:250,:)).^2,1);
    
    % sum of the powers of the two polarizations
    mean_seg_pwr  = mean_seg_x + mean_seg_y ;
    
elseif NumOfChan  < 4
    
    mean_seg_pwr  = mean_seg_x;
    
end

mean_seg_dB = 10*log10(mean_seg_pwr);


%________point_on_TH___________________
figure(101); plot(mean_seg_dB); grid on;
title('Point on lower TH for peak detection','fontsize',18)
xlabel('Samples','fontsize',16)
ylabel('dB','fontsize',16); xlim([1, length(mean_seg_dB)])
clear temp; temp = ginput(1);
TH = temp(2);
close(101)
%______________________________________
if (195<T_RoundTrip && T_RoundTrip<=196)
    N_peak=31;
elseif (490<=T_RoundTrip && T_RoundTrip<=491)
    N_peak=21;
end

done = 'No';
while strcmp(done,'No')
    findpeaks(mean_seg_dB,'MinPeakDistance',min_dis_samp,'NPeaks',N_peak,'MinPeakHeight',TH);
    [a,found_peaks] = findpeaks(mean_seg_dB,'MinPeakDistance',min_dis_samp,'NPeaks',N_peak,'MinPeakHeight',TH);
    title(['N Peaks was set to ', num2str(N_peak), newline, num2str(length(found_peaks)),' were found'],'FontSize',16)
    done = questdlg('\fontsize{13}Correct number of peaks where found?','Peak search', ...
        'Yes','No',struct('Default','Yes','Interpreter','tex'));
    
    if strcmp(done,'No')
        
        ans_modify = questdlg('\fontsize{15}Want to change peaks count OR min distance between peaks?','Peak search', ...
            'Peaks_Count','Min_Distance','Add_manually',struct('Default','Peaks_Count','Interpreter','tex'));
        
        if strcmp(ans_modify,'Peaks_Count')
            prompt = {'Enter Number Of Peaks In The Image'};
            definput = {num2str(N_peak)}; dims = [1 44]; opts.Resize = 'on';
            aa = inputdlg2(prompt,'Peak search',dims,definput,opts);
            N_peak = str2num(aa{1}); % in meters
        elseif strcmp(ans_modify,'Min_Distance')
            prompt = {'Enter Min Distance Between Peaks In Samples'};
            definput = {num2str(min_dis_samp-1)}; dims = [1 44]; opts.Resize = 'on';
            aa = inputdlg2(prompt,'Peak search',dims,definput,opts);
            min_dis_samp = str2num(aa{1}); % in meters
        elseif strcmp(ans_modify,'Add_manually')           
            [aa,~]=ginput(1);           
            [new_Peak_height, hor] = max(mean_seg_dB(round(aa)-10:round(aa)+10));
            found_peaks = [found_peaks , hor + round(aa)-10 - 1];
            a = [a , new_Peak_height];
            done = 'Yes';
        end
        
    end
end
found_peaks = sort(found_peaks,'ascend');
close(gcf)
clear aa ans_modify done a;
a = mean_seg_dB(found_peaks);

figure(101); plot(mean_seg_dB - max(mean_seg_dB)); grid on;
hold on; plot(found_peaks,a- max(mean_seg_dB),'*');
title(['Average power of all segments [dB]',newline,num2str(length(found_peaks)),' peaks'...
    ,newline,'NumOfCopies = ',num2str(NumOfCopies)],'fontsize',18)
xlabel('samples','fontsize',16)
ylabel('Relative power [dB]','fontsize',16);xlim([1, length(mean_seg_dB)])

figure(1001); plot(t_fast*2e8,mean_seg_dB - max(mean_seg_dB)); grid on;
hold on; plot(t_fast(found_peaks)*2e8,a- max(mean_seg_dB),'*');
ax =gca;    ax.XLabel.Interpreter = 'LaTex';
title(['Average power of all segments [dB]',newline,num2str(length(found_peaks)),' peaks'...
    ,newline,'NumOfCopies = ',num2str(NumOfCopies)],'fontsize',18)
str = '$\tilde{z}$ = $v$ ( $\tau$ mod $\tau_{scan}$ ) [m]';
xlabel(str,'fontsize',16)
ylabel('Relative power [dB]','fontsize',16); xlim([1, t_fast(end)*2e8])

figure(2002); plot(t_fast*2e8,mean_seg_pwr); grid on;
ax =gca;    ax.XLabel.Interpreter = 'LaTex';
title(['Average power of all segments [dB]',newline,num2str(length(found_peaks)),' peaks'...
    ,newline,'NumOfCopies = ',num2str(NumOfCopies)],'fontsize',18)
str = '$\tilde{z}$ = $v$ ( $\tau$ mod $\tau_{scan}$ ) [m]';
xlabel(str,'fontsize',16);
ylabel('Relative Power [Optical Watt]','fontsize',16); xlim([1, t_fast(end)*2e8])

%% save m file & figures of power vs. mod(z,L)

fig_name = [num2str(code_length),'b_',num2str(pulse_dur),'ns'];
string_idx = find(fig_name=='.');
fig_name(string_idx(1:end))='_';

if NumOfCopies==1
    savefig(101 ,[results_folder, '\mean_seg_pwr_',fig_name,'_samples'],'compact')
    savefig(1001,[results_folder, '\mean_seg_pwr_',fig_name,'_meter'],'compact')
    savefig(2002,[results_folder, '\mean_seg_pwr_',fig_name,'_meter_lin'],'compact')
    close (2002);
else
    savefig(101 ,[results_folder, '\mean_seg_pwr_',fig_name,'_',num2str(NumOfCopies),'_copies'],'compact')
    savefig(1001 ,[results_folder, '\mean_seg_pwr_',fig_name,'_',num2str(NumOfCopies),'_meter'],'compact')
end
close (1001); 

%% find peaks per polarization
if NumOfChan == 4
    
    for ii=1:length(found_peaks)
        pk_zone = found_peaks(ii)-3 : found_peaks(ii)+3 ;
        if pk_zone(end)>length(mean_seg_x)
            pk_zone = found_peaks(ii)-3 : length(mean_seg_x);
        end
        [ ~ , idx_x]  = max(mean_seg_x(pk_zone));
        [ ~ , idx_y]  = max(mean_seg_y(pk_zone));
        pk_x(ii) = pk_zone(idx_x);
        pk_y(ii) = pk_zone(idx_y);
    end
    
    % plot the 2 pol on the same graph
    mean_seg_x_db = 10*log10(mean_seg_x);
    mean_seg_y_db = 10*log10(mean_seg_y);
    
    max_val =  max(max(mean_seg_x_db(pk_x)),max(mean_seg_y_db(pk_y))) ;
    figure(11);
    plot(mean_seg_x_db-max_val);hold on; plot(pk_x,mean_seg_x_db(pk_x)-max_val,'r*'); grid on;
    hold on; plot(mean_seg_y_db-max_val);hold on; plot(pk_y,mean_seg_y_db(pk_y)-max_val,'g*');
    title(['average power for the 2 polarizations [dB]',newline,num2str(length(pk_y)),' peaks'],'fontsize',18)
    xlabel('samples','fontsize',16)
    ylabel('dB','fontsize',16);xlim([1, length(mean_seg_dB)])
    legend('x pol','','y pol','')
    ylim([min(mean_seg_x_db)-max_val , 0.1]);
    savefig(11, [results_folder, '\mean_seg_pwr_of_2_polarizations'],'compact');
       
    
    for ii=1:length(found_peaks)
        pk_zone = found_peaks(ii)-3 : found_peaks(ii)+3 ;
        if pk_zone(end)>length(mean_seg_x)
            pk_zone = found_peaks(ii)-3 : length(mean_seg_x);
        end
        [ ~ , idx_x]  = max(mean_seg_x(pk_zone));
        [ ~ , idx_y]  = max(mean_seg_y(pk_zone));
        pk_x(ii) = pk_zone(idx_x);
        pk_y(ii) = pk_zone(idx_y);
    end
    
end


%% mark pulse width (in samples) so we can remove the pulses out and left with the noise floor (29/09/22) :

window_size=[];
% arbitrary pick a pulse:
peak_idx = randi(length(found_peaks));
% open display over a single pulse:
Half_peak_width = round((length(mean_seg_dB)/length(found_peaks))/2);
if (found_peaks(peak_idx)-Half_peak_width) < 0
    Half_peak_width = found_peaks(peak_idx)-1;
end
if (found_peaks(peak_idx)+Half_peak_width > length(mean_seg_dB))
    while (isempty(window_size))
        % plot pulse and mark its total width
        figure(); plot(mean_seg_dB( found_peaks(peak_idx)-Half_peak_width  : end )); grid on;
        xlabel('samples','FontSize',16); title(['Noise floor calc',newline,'click twice, one on each side of the pulse'...
            ,newline, ' the area between the two clicks will be removed for noise floor calc'])
        legend(['Peak # ',num2str(peak_idx)])
        UserInput = ginput(2);
        Left_samp =  round(min(UserInput(1,1), UserInput(2,1)));
        Right_samp = round(max(UserInput(1,1), UserInput(2,1)));
        hold on; plot(Left_samp: Right_samp , mean_seg_dB( found_peaks(peak_idx) - Half_peak_width + Left_samp -1   : found_peaks(peak_idx) - Half_peak_width + Right_samp -1  ),'r*');
        PulseWidth= Right_samp - Left_samp;  % [Samples]
        text(round(Left_samp+PulseWidth/2), mean_seg_dB(found_peaks(peak_idx))+1,['Pulse Width is ',num2str(PulseWidth) ,' Samples'],'FontSize',16 )
        
        % set the desired pulse width to remove:
        prompt = {['Enter window size for peak removal:',newline,'cancle will bring back the interactive figure']};
        definput = {num2str(PulseWidth)}; dims = [1 40]; opts.Resize = 'on';
        window_size = inputdlg2(prompt,'Legendre bits',dims,definput,opts);
        if ~isempty(window_size)
            window_size = str2num(window_size{1}); % num of samples to discard per peak (for noise calc)
        else
            close(gcf)
        end
    end
       
else

    while (isempty(window_size))
        % plot pulse and mark its total width
        figure(); plot(mean_seg_dB( found_peaks(peak_idx)-Half_peak_width  :  found_peaks(peak_idx)+Half_peak_width )); grid on;
        xlabel('samples','FontSize',16); title(['Noise floor calc',newline,'click twice, one on each side of the pulse'...
            ,newline, ' the area between the two clicks will be removed for noise floor calc'])
        legend(['Peak # ',num2str(peak_idx)])
        UserInput = ginput(2);
        Left_samp =  round(min(UserInput(1,1), UserInput(2,1)));
        Right_samp = round(max(UserInput(1,1), UserInput(2,1)));
        hold on; plot(Left_samp: Right_samp , mean_seg_dB( found_peaks(peak_idx) - Half_peak_width + Left_samp -1   : found_peaks(peak_idx) - Half_peak_width + Right_samp -1  ),'r*');
        PulseWidth= Right_samp - Left_samp;  % [Samples]
        text(round(Left_samp+PulseWidth/2), mean_seg_dB(found_peaks(peak_idx))+1,['Pulse Width is ',num2str(PulseWidth) ,' Samples'],'FontSize',16 )
        
        % set the desired pulse width to remove:
        prompt = {['Enter window size for peak removal:',newline,'cancel will bring back the interactive figure']};
        definput = {num2str(PulseWidth)}; dims = [1 40]; opts.Resize = 'on';
        window_size = inputdlg(prompt,'Legendre bits',1,definput,opts.Resize);

        if ~isempty(window_size)
            window_size = str2num(window_size{1}); % num of samples to discard per peak (for noise calc)
        else
            close(gcf)
        end
    end
end
close(gcf)
%% calc static SNR
[snrx,snry] = max_polarization(figure(11),window_size,scan_num);
[combined_snr] = combined_max_polarization(figure(101),window_size,scan_num);

save([file(1).folder,'\SNR'],'snrx','snry','combined_snr')

%% measure peak temporal width
if (pk_x(1)-100>0)
    jj=1;
else
    jj=2;
end
single_peak = mean_seg_x(pk_x(jj)-space  :pk_x(jj)+space);

aa = find(single_peak >= (max(single_peak)/2));
peak_width_samp = length(aa);

if debug
    figure; plot(single_peak,'linewidth',2); grid on;
    hold on; plot([aa(1),aa(end)], [single_peak(aa(1)),single_peak(aa(end))],'*')
    title(['Peak full width at half max is ', num2str(peak_width_samp),' samples',newline,'~ ',num2str(peak_width_samp*dt*1e9),' ns']);
    xlabel('samples','FontSize',16)
end
%% save original mat for reference
c_data_x_orig = c_data_x;

%% rearrange the waterfall map
if 0
    if T_RoundTrip > 490 % only for a decending order array
        
        % move the reflections in the slow axis so that each line will be the
        % result of the same interrogating pulse
        % figure; imagesc(abs(c_data_x_orig))
        T_scan = code_length*pulse_dur; % [ns]
        if (T_scan > T_RoundTrip) && (T_scan < 2*T_RoundTrip)
            
            c_data_x_arranged = zeros(size(c_data_x));
            c_data_y_arranged = zeros(size(c_data_y));
            
            for kk=1:length(found_peaks)
                temp =  circshift( c_data_x_orig(:,found_peaks(kk)-ceil(peak_width_samp/2):found_peaks(kk)+ceil(peak_width_samp/2)) , -(length(found_peaks)-kk))  ;
                c_data_x_arranged(:,found_peaks(kk)-ceil(peak_width_samp/2):found_peaks(kk)+ceil(peak_width_samp/2)) = temp;
            end
            c_data_x_arranged(end-length(found_peaks)+1:end,:)=0;
            
            if NumOfChan==4
                for kk=1:length(found_peaks)
                    temp =  circshift( c_data_y(:,found_peaks(kk)-ceil(peak_width_samp/2):found_peaks(kk)+ceil(peak_width_samp/2)) , -(length(found_peaks)-kk))  ;
                    c_data_y_arranged(:,found_peaks(kk)-ceil(peak_width_samp/2):found_peaks(kk)+ceil(peak_width_samp/2)) = temp;
                end
                c_data_y_arranged(end-length(found_peaks)+1:end,:)=0;
            end
            
            figure; imagesc(abs(c_data_x_arranged));
            
            c_data_x = c_data_x_arranged;
            c_data_y = c_data_y_arranged;
            
            clear c_data_x_arranged c_data_y_arranged temp
        end
        
    else
        c_data_x_arranged = zeros(size(c_data_x));
        c_data_y_arranged = zeros(size(c_data_y));
        
        NumPeaks = length(found_peaks)-2;
        for kk=1+2:length(found_peaks)
            temp =  circshift( c_data_x_orig(:,found_peaks(kk)-ceil(peak_width_samp/2):found_peaks(kk)+ceil(peak_width_samp/2)) , +(NumPeaks-kk))  ;
            c_data_x_arranged(:,found_peaks(kk)-ceil(peak_width_samp/2):found_peaks(kk)+ceil(peak_width_samp/2)) = temp;
        end
        c_data_x_arranged(end-NumPeaks+1:end,:)=0;
        c_data_x = c_data_x_arranged;
        
    end
end
%% plot abs value of complex data in time
if debug
    figure(99)
    imagesc (t_fast*2e8 , t_slow*1e3, abs(c_data_x) - max(abs(c_data_x(:))));
    colorbar
    ylabel('slow time [ms]','fontsize',16)
    xlabel('Relative distance [m]','fontsize',16)
    title('abs of detected data [dB] ','fontsize',18)
    
    saveas(99,[results_folder, '\waterfall_pwr'],'png')
    close(99);
end
%% design HP filter that will work on the slow axis
if debug
    f_pass = 50e3; % Hz
    f_pass = f_pass/(Fs_slow/2);
    A_stop = 50;
    A_pass = 1;
    filter_order = size(c_data_x);
    filter_order = filter_order(1);
    filter_order = round(filter_order/4);
    
    d = fdesign.highpass('N,Fp,Ast,Ap',filter_order,...
        f_pass,A_stop,A_pass);
    Hd = design(d,'equiripple'); % Design Methods:  ellip / equiripple { ifir / kaiserwin / butter / cheby1 }
    filter_coeff = Hd.Numerator;
    
    if (debug)
        [h,f] = freqz(Hd,[],Fs_slow);
        figure;subplot(2,1,1)
        plot(f, dB20(h)); grid
        xlabel('freq [Hz]'); ylabel('20*log_1_0');
        title('highpass filter','fontsize',14)
        
        subplot(2,1,2); plot(filter_coeff);
        xlabel('samples in the slow axis (or scan number)')
        title('filter coeff','fontsize',14); grid;
    end
    
end

%% extract diff phase between two consecutive reflectors

concat_spools=0;

if concat_spools==1
    % change peaks order according to their real position order in the lab
    % for the combination: 2 -> acoustic bath ->16
    pairs_order  = [2,1,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3];
    msgbox(['this code is for a spesific peaks order',newline,'check "pairs_order" variable'])
elseif concat_spools==0
    % for a single spool interrogation
    pairs_order  = [length(found_peaks):-1:1];
end

snrx_oreder = snrx(pairs_order);
snry_oreder = snry(pairs_order);
found_peaks_order = found_peaks(pairs_order);

[best_pol,diff_phase] = extract_diff_phase(SNRTH,snrx_oreder,snry_oreder,found_peaks_order,c_data_x,c_data_y,Fs_slow,results_folder,0);

%% delete HEAVY vector/matrice

clearvars I_data_x I_data_y c_data_x_orig c_data_db_x Q_data_x Q_data_y

%% detected (differential) phases in time domain

if debug
    figure;
    plot(t_slow*1e3, diff_phase*1e3);
    xlabel('Time [ms]'); ylabel('Phase [mrad]'); title('diff phase between consecutive pulses')
    legend; grid on;
    savefig(gcf,[results_folder, '\diff_phase_time_domain'],'compact')
    close;
end
%% PSD of diffrential phase

T_window = 0.9; % ms
prompt = {['Insert window size for PWELCH in [ms]:']};
definput = {num2str(T_window)}; dims = [1 40]; opts.Resize = 'on';
T_window = str2num(cell2mat(inputdlg2(prompt,'PSD',dims,definput,opts))); % ms

window_size = round(T_window*1e-3*Fs_slow); % samples
[pxx , f] = pwelch (diff_phase,window_size,window_size/10,[],Fs_slow);

f_10_idx = find ( f > 10e3,1,'first');  % index of f = 10kHz

figure(83);
plot( f/1e3 , sqrt(pxx)*1e3)
ylabel('mrad / sqrt(Hz)'); xlabel('frequency [kHz]');
ylim([0 , max(max(sqrt(pxx(f_10_idx:end,:))*1e3))])

savefig(83,[results_folder, '\PSD_linear'],'compact')
close (83)

figure(183);
plot( f/1e3 , 10*log10((pxx)))
ylabel('10*log_1_0 rad / sqrt(Hz)'); xlabel('frequency [kHz]');

savefig(183,[results_folder, '\PSD_log_scale'],'compact')
close (183)

%% design BP filter that will work on the slow axis
aa = questdlg('\fontsize{13}Perform BPF on phase data? ','BandPass', ...
        'Yes','No',struct('Default','Yes','Interpreter','tex'));
            
if strcmp(aa,'Yes')
    
    temp = find (file(1).name == 'k');
    temp = temp(end);
    f_acoustic_char = file(1).name(temp-3:temp-1);
    
    if f_acoustic_char(1)=='_'
        f_acoustic_char(1)='';
    end
    
    prompt = {'Enter acoustic frequency [kHz]'};
    definput = {f_acoustic_char}; dims = [1 44]; opts.Resize = 'on';
    aa = inputdlg2(prompt,'underwater acoustics',dims,definput,opts);
    f_acoustic = str2num(aa{1})*1e3;
    clear aa
    
    f_stop_1 = f_acoustic-10*1e3;
    f_pass_1 = f_acoustic-5*1e3;
    f_pass_2 = f_acoustic+5*1e3;
    f_stop_2 = f_acoustic+10*1e3;
    
    A_stop_1 = 50;
    A_pass = 2;
    A_stop_2 = 50;
    
    d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',f_stop_1,...
        f_pass_1 , f_pass_2,f_stop_2,A_stop_1,A_pass, A_stop_2, Fs_slow);
    Hd = design(d,'butter'); % Design Methods:  ellip / equiripple { ifir / kaiserwin / butter / cheby1 }
    
    %     fvtool(Hd) % filter_coeff = Hd.Numerator;
    
    diff_phase_filtered = filter(Hd , diff_phase(:,1:end));

    figure(2021); plot(t_slow*1e3,diff_phase_filtered(:,15:end)*1e3); legend;
    grid on; ylabel('phase [mrad]'); title(['Diff phase between consecutive pulses', newline , 'after BPF']);
    xlabel('Time [ms]')
    savefig(2021,[results_folder, '\diff_phase_time_domain_filtered'],'compact')
    close (2021)
    
    figure(2022); imagesc([],t_slow*1e3,abs(diff_phase_filtered(:,15:end)).^2);
    ylabel(['time [ms]']); caxis([0,0.4*max(abs(diff_phase_filtered(:)).^2)]);
    xlabel('Mandrel number')
    savefig(2022,[results_folder, '\diff_phase_time_domain_waterfall'],'compact')
    close (2022)
    
else
    
    diff_phase_filtered = diff_phase;
    
end
clear aa;

%% UnderWater communication using interrogation pulsr

UWC = questdlg('\fontsize{16}Is it a communication exp? ','AutoCorr', ...
    'Yes','No',struct('Default','Yes','Interpreter','tex'));

if strcmp(UWC,'Yes')
    
    % acoustic data was transmitted as follow:
    %_____pulse______15ms__________data_(bpsk)____
    
    aa = questdlg('\fontsize{16}chose bit rate [kbps] ','Bit Rate', ...
        '4','8','other',struct('Default','4','Interpreter','tex'));
    BitRate = str2num(aa); % [kbps]
    
    if isempty(BitRate)
        aa = questdlg('\fontsize{16}chose bit rate [kbps] ','Bit Rate', ...
            '1','2','other',struct('Default','2','Interpreter','tex'));
        BitRate = str2num(aa); % [kbps]
    end
    
    aa = questdlg('\fontsize{16}chose number of transmitted bits','Number of bits', ...
        '8','32','other',struct('Default','4','Interpreter','tex'));
    N_bits = str2num(aa); % [kbps]
    
    for Mandrel_idx = 1:size(diff_phase_filtered,2)
        Mandrel = diff_phase_filtered(:,Mandrel_idx);
        [~,ssb,ff] = my_fft_func(Mandrel,Fs_slow);%     figure; plot(ff, abs(ssb))
        Mandrel_BB = (ifft(ifftshift(ssb))).*exp(-1i*2*pi*f_acoustic.*t_slow.');
        AutoCorr(:,Mandrel_idx) = xcorr(Mandrel_BB,Mandrel_BB);
        
        if Mandrel_idx==1
            laa = length(AutoCorr(:,Mandrel_idx));
            T_slow=laa*dt_slow;
            t_slow1 = -T_slow/2:dt_slow:(T_slow/2 - dt_slow);
            sample_point = 15e-3 + [0:1:N_bits-1]*(1/(BitRate*1e3));
            for ii = 1:length(sample_point)
                idx_4_MF(ii) = find (t_slow1<=sample_point(ii),1,'last');
            end
        end
    end
    
    % save autocorrelation matrix:
    save ([results_folder,'\AutoCorr_',num2str(N_bits),'bit_',num2str(BitRate),'kbps'],'AutoCorr','t_slow1','sample_point','idx_4_MF')
    
    Manrel_idx = 8;
    figure(12); plot(t_slow1 ,real(AutoCorr(:,Manrel_idx)))
    grid on;
    hold on; plot(t_slow1(idx_4_MF) ,real(AutoCorr(idx_4_MF,Manrel_idx)),'o','linewidth',3.6)
    for kk = 1:N_bits
        figure(12); hold on;
        line([sample_point(kk) sample_point(kk)],[min(real(AutoCorr(:,Manrel_idx))) max(real(AutoCorr(:,Manrel_idx)))]);
    end
    hold on; line([sample_point(1) sample_point(end)],[0 0], 'LineWidth',3 ,'Color','black');
    savefig(12,[results_folder, '\AutoCorr_Mandrel_',num2str(Manrel_idx)],'compact')
    
else
    
    %% Impulse response or Matched Filtering
    
    % % Down conversion to BB:
    % Mandrel = diff_phase_filtered(:,10);
    % [~,ssb,ff] = my_fft_func(Mandrel,Fs_slow);%     figure; plot(ff, abs(ssb))
    % Mandrel_BB = (ifft(ifftshift(ssb))).*exp(-1i*2*pi*f_acoustic.*t_slow.');
    
    aa = questdlg('Perform Matched filter on phase data? ','Matched Filter', ...
        'Yes','No','Yes');
    
    if strcmp(aa,'Yes')
        
        % pick a mandrel:
        Mandrel_idx = 10;
        
        % loda MF:
        load(['D:\Nadav_D\experiments\FBG_array_2022\FBG_20m\0km_spool\communication\first_try\Impulse_Response\2\MatchFilter_lst'])
        MF2 = MF_BB(:,Mandrel_idx);
        
        load(['D:\Nadav_D\experiments\FBG_array_2022\FBG_20m\0km_spool\communication\first_try\Impulse_Response\4\MatchFilter_lst'])
        MF4 = MF_BB(:,Mandrel_idx);
        
        % Down conversion to BB:
        Mandrel = diff_phase_filtered(:,Mandrel_idx);
        [~,ssb,ff] = my_fft_func(Mandrel,Fs_slow);%     figure; plot(ff, abs(ssb))
        Mandrel_BB = (ifft(ifftshift(ssb))).*exp(-1i*2*pi*f_acoustic.*t_slow.');
        
        MF_out = xcorr(Mandrel_BB ,MF4 );
        
        figure(88); plot([-(length(MF_out)/2):1:(length(MF_out)/2-1)]*(1/Fs_slow)*1e3 , real(MF_out(1:end)))
        MF_out = xcorr(Mandrel_BB ,MF2 );
        hold on;
        plot([-(length(MF_out)/2):1:(length(MF_out)/2-1)]*(1/Fs_slow)*1e3 , real(MF_out(1:end)))
        grid on; xlabel('TIME [msec]'); ylabel(['lin']); title('After MF')
        savefig(88,[results_folder, '\MF_result'],'compact')
        
        save([file(1).folder,'\diff_phase_filtered'],'diff_phase_filtered');
        
    elseif  strcmp(aa,'No')
        aa = questdlg('Save data as Matched Filter? ','Matched Filter', ...
            'Yes','No','Yes');
        
        if strcmp(aa,'Yes')
            % move Matched Filter to BB
            MF_lst = diff_phase_filtered;
            
            [a,ssb,ff] = my_fft_mat_func(MF_lst,Fs_slow); close;
            %         figure; plot(ff, abs(ssb(:,end-1)))
            
            MF_BB = (ifft(ifftshift(ssb))).*exp(-1i*2*pi*f_acoustic.*t_slow.');
            
            save([file(1).folder,'\MatchFilter_lst'],'MF_lst','MF_BB');
            clear MF_lst;
        end
        save([file(1).folder,'\diff_phase_filtered'],'diff_phase_filtered');
    end
    
end    
    
    %% Spectrogram
    aa = questdlg('Do u want to plot spectrogram? ','Spectrogram', ...
        'Yes','No','Yes');
    aa = 'No';
    if strcmp(aa,'Yes')
        
        SegmentUnderTest=diff_phase_filtered(:,end);
        Segment_before_last = diff_phase_filtered(:,end-1);
        SegmentInDrum = diff_phase_filtered(:,end-2);
        
        figure; plot(t_slow*1e3 ,SegmentUnderTest);
        hold on; plot(t_slow*1e3 ,Segment_before_last);
        hold on; plot(t_slow*1e3 ,SegmentInDrum);
        xlabel('Time  [ms]'); ylabel('Phase [rad]');legend('last','before last' ,'inside drum')
        
        prompt = {'Insert temporal window duration in ms'};
        definput = {'0.5'}; dims = [1 38]; opts.Resize = 'on';
        aa = inputdlg2(prompt,'Spectrogram',dims,definput,opts);
        Win_duration = str2num(aa{1})*1e-3; % in sec
        clear aa
        
        WindowSize = round(Win_duration*Fs_slow);
        noverlap = round(WindowSize*0.95);
        nfft = WindowSize/2;
        
        [s,f,t] = spectrogram(SegmentUnderTest, WindowSize, noverlap,[], Fs_slow);
        figure(28);imagesc(f/1e3 , t*1e3 ,10*log10(abs(s.')))
        xlabel('Freq [kHz]'); ylabel('Time [ms]'); colorbar;
        title('Last segment'); ax = gca; ax.FontSize=16;
        savefig(28,[results_folder, '\spectrogram_last_segment'],'compact')
        
        % find max freq per row and plot it:
        %     max_per_row = max(s);
        %     s_nrom = s./max_per_row;
        %     figure(82);imagesc(f/1e3 , t*1e3 ,10*log10(abs(s_nrom.')))
        %     xlabel('Freq [kHz]'); ylabel('Time [ms]'); colorbar;
        %     title('Last segment'); ax = gca; ax.FontSize=16;
        
        
        [s,f,t] = spectrogram(Segment_before_last, WindowSize, noverlap,[], Fs_slow);
        figure(29);imagesc(f/1e3 , t*1e3 ,10*log10(abs(s.'))); xlim([0 , 100])
        xlabel('Freq [kHz]'); ylabel('Time [ms]'); colorbar
        title('Before last segment'); ax = gca; ax.FontSize=16;
        savefig(29,[results_folder, '\spectrogram_before_last_segment'],'compact')
        close (29)
        
        [s,f,t] = spectrogram(diff_phase_filtered(:,end-4), WindowSize, noverlap,[], Fs_slow);
        figure(29);imagesc(f/1e3 , t*1e3 ,10*log10(abs(s.'))); xlim([5 , 80])
        xlabel('Freq [kHz]'); ylabel('Time [ms]'); colorbar
        title('Before last segment'); ax = gca; ax.FontSize=16;
        %     savefig(30,[results_folder, '\spectrogram_before_last_segment'],'compact')
        %     close (30)
        
        figure; plot(t_slow*1e3 ,diff_phase_filtered(:,end-4));
        xlabel('Time  [ms]'); ylabel('Phase [rad]');
        
    end
    %% peaks amplitude analysis
    if debug
        for N=1:length(found_peaks_order)
            
            amp_data_x = abs(c_data_x(:,found_peaks_order(N)));
            
            %     figure;plot(t_slow*1e6,amp_data_x, 'linewidth',2);
            %     grid on; xlabel('time [us]')
            
            [~,ssb_amp,ff] = my_fft_func(amp_data_x,Fs_slow);
            
            if N==1
                figure(105);
            end
            
            plot(ff/1e3,20*log10(abs(ssb_amp)),'-','linewidth',1.5);grid on
            hold on;
        end
        
        xlim([0,ff(end)/NumOfCopies]/1e3);
        title(['detected peaks amplitude spectrum ',newline,'SSB',newline,'X polarization'],'fontsize',18)
        xlabel('freq [kHz]','fontsize',16)
        ylabel('dB','fontsize',16);
        legend;
        savefig(105,[results_folder, '\amp_freq_domain_all_peaks_Xpol'],'compact')
        close (105)
    end
    %% done
    disp('Done');
    close all;