% written by Nadav A - March 14
% modified: November 2022, Yair BN
% the code generate a Legendre and trernary sequence 
% and plots them

clc; close all; clear;
%debug = 1;
%% _____user_params_____________________________
code_length = 59; % must be a prime number 
while ~isprime(code_length) || mod(code_length,4) ~= 3
    warning(['invalid input', newline, 'Please enter a prime number'])
    prompt = "Please enter the code length (should be prime): ";
    code_length = input(prompt);
end

%_____pulse_params_____________________________
Fs = 10;  
Fs = Fs*1e9;
dt =  1/Fs;

pulse_dur = 3.5;                 %[nsec]
T_code = pulse_dur*code_length; %[nsec]

pulse_duration = ( pulse_dur*1e-9 );  %[sec]
%_____time_vector________________________________
t = -pulse_duration/2 : dt : pulse_duration/2-dt;
%_____pulse_shape_____________________________
pulse = cos(2*pi*t/pulse_duration) + 1 ;
pulse = [zeros(1,2*16-length(pulse)),pulse];

%_____data_____________________________
data = imag(perfect_periodic_Legendre_waveform(code_length));
data = (data/(max(data)))*2-1;
close

waveform_leg = kron(data,pulse)/2;
time_vector_leg = 0:dt:(length(waveform_leg)-1)*dt;
data_p = rectpulse(data,length(pulse))*max(waveform_leg);

%% ___________Ternary___________________________

% load the desired ternary code

Ter_data = [1, 1, 1, 1, 1, -1, 1, 0, -1, 1, -1, 0, -1, -1, -1, 1, -1, -1, 0, 1, -1, 1, 1, -1, -1, 1, 1, 0, -1, 1, -1, 1, 0, 0, 1, 0, 1, 1, 1, 1, -1, -1, -1, 1, 1, 0, 1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1];
code_length = length(Ter_data);

%% _____pulse_params_____________________________

Fs = 10;  
Fs = Fs*1e9;
dt =  1/Fs;

pulse_dur = 20;                                          %[nsec]
pulse_dur = 3.5; %[nsec]
T_code = pulse_dur*code_length;                          %[nsec]

pulse_duration = ( pulse_dur*1e-9 );  %[sec]
%_____time_vector________________________________
t = -pulse_duration/2 : dt : pulse_duration/2-dt;
%_____pulse_shape_____________________________
pulse = cos(2*pi*t/pulse_duration) + 1 ;
pulse = [zeros(1,2*16-length(pulse)),pulse];

waveform_ter = kron(Ter_data,pulse)/2;
time_vector_ter = 0:dt:(length(waveform_ter)-1)*dt;

% data_p = rectpulse(Ter_data,length(pulse))*max(waveform);

figure(9);

subplot(2,1,1);
plot(time_vector_ter*1e9 , waveform_ter,'r-' ); hold on;
title('Ternary code','fontsize',26);
xlabel('Time [ns]','fontsize',26,'FontWeight','bold');
ylabel('E _{Out}','fontsize',26,'FontWeight','bold');
xlim([21 ,56]);
grid on; box off;
ax = gca;
ax.FontSize = 28
ax.Children.LineWidth = 5
yline(0,'--k','LineWidth',3)

subplot(2,1,2);
plot(time_vector_leg*1e9 , waveform_leg,'r-' );
title('Legendre Code','fontsize',26);
xlabel('Time [ns]','fontsize',26,'FontWeight','bold')
ylabel('E _{Out}','fontsize',26,'FontWeight','bold')
xlim([21,56])
grid on; box off;
ax = gca;
ax.FontSize = 28
ax.Children.LineWidth = 5
yline(0,'--k','LineWidth',3)

%% __________save_to_file_______________________
fig_name = 'Oasis_Ter_and_leg_plot.fig';
savefig(9 ,['D:\alpha\Yair\MATLAB\Figures\',fig_name])

%% plot the last 10 bits
tau_ns_ter = 14; %[ns] 
tau_s_ter = tau_ns_ter*1e-9; %[sec] 

tau_ns_leg = 21; %[ns] 
tau_s_leg = tau_ns_leg*1e-9; %[sec] 

figure(10);

subplot(2,1,1);
plot(time_vector_ter(end-tau_s_ter*Fs:end)*1e9 , waveform_ter(end-tau_s_ter*Fs:end),'r-' ); hold on;
title('Ternary code','fontsize',26);
xlabel('Time [ns]','fontsize',26,'FontWeight','bold');
xlim([174 ,199.5]); ylim([-1,+1])
grid on; box off;
ax = gca;
ax.FontSize = 28
ax.Children.LineWidth = 5
y = zeros(1, length(time_vector_leg(end-tau_s_leg*Fs:end)*1e9 )-25);
line(time_vector_leg(end-tau_s_leg*Fs)*1e9:time_vector_leg(end-tau_s_leg*Fs)*1e9+length(y)-1,y,...
    'LineWidth',3,'Color','black','LineStyle','--')
set(ax , 'YColor','none')  % set(ax ,'ytick',[])
aa = str2num(cell2mat([ax.XTickLabel ; {'199'}]));
ax.XTick = [aa];

subplot(2,1,2);
plot(time_vector_leg(end-tau_s_leg*Fs:end)*1e9 , waveform_leg(end-tau_s_leg*Fs:end),'r-' );
title('Legendre Code','fontsize',26);
xlabel('Time [ns]','fontsize',26,'FontWeight','bold')
xlim([174 ,199.5]); ylim([-1,+1])
grid on; box off;
ax = gca;
ax.FontSize = 28;
ax.Children.LineWidth = 5;
% yline(0,'--k','LineWidth',3)
y = zeros(1, length(time_vector_leg(end-tau_s_leg*Fs:end)*1e9 )-25);
line(time_vector_leg(end-tau_s_leg*Fs)*1e9:time_vector_leg(end-tau_s_leg*Fs)*1e9+length(y)-1,y,...
    'LineWidth',3,'Color','black','LineStyle','--')
set(ax , 'YColor','none')
aa = str2num(cell2mat([ax.XTickLabel ; {'199'}]));
ax.XTick = [aa];% ax.XTickLabel


%% __________save_to_file_______________________
fig_name = 'Oasis_Ter_and_leg_plot_end.fig';
savefig(10 ,['D:\alpha\Yair\MATLAB\Figures\',fig_name])
%close all;


