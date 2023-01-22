function s = perfect_periodic_Legendre_waveform(N)
% Generates a periodic coded signal using 2 or 3 phases.
% The signal exhibits perfect peiodic autocorelation 
% N is any odd prime
Nspt=sprintf('%g element phase-coded waveform ', N);
s=ones(1,N);
if isprime(N)==0
    disp('Not a prime')
    return
end
if rem((N+3)/4, 1)==0
    c = 0.25*(N-1);
    c1 = 2*1/c-1/(2*c^2);
    c2 = 2*1/c-1/(4*c^2);
    arg2=acos(-c1/2-sqrt((c1/2)^2-c2));
    s(mod((1:N-1).^2,N)+1) = exp(1i*arg2);
    s(1) = exp(1i*arg2/2);
else
    arg3=acos(-(N-1)/(N+1));
    s(mod((1:N-1).^2,N)+1) = exp(1i*arg3);
end
d = abs(ifft(fft(s).*conj(fft(s))));
figure; plot(d, 'k')
title(['Periodic autocorrelation of ', Nspt])
end

    
    
    
 
