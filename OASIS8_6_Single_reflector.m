close all; clc; clearvars
%%
%% load legendre codes
leg_1_length = input('Enter first legendre code length:')

leg_1_data = imag(perfect_periodic_Legendre_waveform(leg_1_length));
leg_1_data = (leg_1_data/(max(leg_1_data)))*2-1;
%% simulation params
Fs = 20e9;
dt =  1/Fs;

%% _____first_leg_pulse_params_____________________________


pulse_dur = input('Insert first Legendre pulse duration in nano-sec:')              %[nsec]
T_code_leg_1 = pulse_dur*leg_1_length; %[nsec]

pulse_duration = ( pulse_dur*1e-9 );  %[sec]
%_____time_vector________________________________
t = -pulse_duration/2 : dt : pulse_duration/2-dt;
%_____pulse_shape_____________________________
pulse = cos(2*pi*t/pulse_duration) + 1 ;

waveform_leg_1 = kron(leg_1_data,pulse);

L_1_AR = abs(xcorr(waveform_leg_1,[waveform_leg_1,waveform_leg_1,waveform_leg_1],2*length(waveform_leg_1)));
L_1_AR = circshift(L_1_AR(1:2*length(waveform_leg_1)),-260);
%%
close all;
Folder = 'D:\alpha\Yair\OASIS_GRAPHS\Reflection\';
h1 = openfig([Folder , 'Leg59_PeriodicAutoCorr.fig']);

hold on; plot([1:length(L_1_AR)]*40/length(L_1_AR),20*log10((L_1_AR)/max(L_1_AR)),'Linewidth',2)
ylim([-52,0.5])
ax =gca;
ax.FontSize=25;
ax.Title.FontSize=32;
Lines = h1.Children.Children
legend('Measured','Simulation');

Lines(1).Color=[0.8,0.1,0.1];
Lines(2).Color=[0.1,0.1,0.1];
ax.XLim=[1,36];

uistack(Lines(2),'top')
%% save figure

savefig(gcf,[Folder, 'Leg59_PeriodicAutoCorr_edit.fig'])