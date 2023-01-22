% written by Nadav A - March 14
% modified: April 2022
% the code generate a Legendre sequence 
% write it into txt file
% this file is then inserted to the AWG

clc; close all; clear;
debug = 1;
%% _____user_params_____________________________
code_length = 59; % must be a prime number 
while ~isprime(code_length) || mod(code_length,4) ~= 3
    warning(['invalid input', newline, 'Please enter a prime number'])
    prompt = "Please enter the code length (should be prime): ";
    code_length = input(prompt);
end

%_____pulse_params_____________________________


Fs = input('Enter sampling frequency (In units of GHz): ');  
Fs = Fs*1e9;
dt =  1/Fs;

pulse_dur = input("Enter Pulse Duration:");                 %[nsec]
T_code = pulse_dur*code_length; %[nsec]

pulse_duration = ( pulse_dur*1e-9 );  %[sec]
%_____time_vector________________________________
t = -pulse_duration/2 : dt : pulse_duration/2-dt;
%_____pulse_shape_____________________________
pulse = cos(2*pi*t/pulse_duration) + 1 ;
pulse = [zeros(1,2*16-length(pulse)),pulse];

figure(1);
plot(t*1e9,abs(pulse)); hold on ; plot(t*1e9,real(pulse),'*'); legend('abs','real');
xlabel('time [ns]','fontsize',12);
xlim([t(1)-dt,t(end)+dt]*1e9); grid
%_____data_____________________________
data = imag(perfect_periodic_Legendre_waveform(code_length));
data = (data/(max(data)))*2-1;
close

figure(2)
stem( (data)); grid on;
title('legendre code','fontsize',16);
xlabel('samples','fontsize',16)
xlim([-0.5 ,  code_length+0.5])
ylim([-1-0.5 ,  1+0.5])

%% _pulse_shping_on_data_________________________
% a code sequence will look :
% |<----L----->|<----L----->|<----L----->|

%%%%%___make_new_pulse_with_new_pulse_duration____%%%%%%%%%
% pulse_dur = T_code/(code_length/DC);   %[nsec]
% pulse_duration = ( pulse_dur*1e-9 );  %[sec]
% %_____time_vector________________________________
% t = -pulse_duration/2 : dt : pulse_duration/2-dt;
% %_____pulse_shape_____________________________
% pulse = cos(2*pi*t/pulse_duration) +1 ;
% 
% %%%%%%%%%%%%%%%%%%%
% temp = zeros ([ padding , size(data,2) ]);
% data2 = [data , temp] ;
% code_length_w_padd = size(data2,2);


waveform = kron(data,pulse);

data_p = rectpulse(data,length(pulse))*max(waveform);

if (debug == 1)
    figure
    plot( waveform,'--o' ); hold on;
    plot( data_p,'--o' );grid on; 
    title('Legendre code','fontsize',16);
    xlabel('samples','fontsize',16)
    xlim([-0.5 ,  code_length*length(pulse)+0.5])
end

%% __________save_to_file_______________________

FileName = ['Leg_',num2str(code_length),'_',num2str(pulse_dur),'ns.txt']
fileID = fopen([cd,'\txt_waves\',FileName],'w');

dlmwrite([cd, '\txt_waves\' ,FileName],waveform.','precision','%2.10f') % fprintf(fileID,'%2.9f',waveform);

fclose(fileID);
    
% my_fft_func(waveform,Fs);


