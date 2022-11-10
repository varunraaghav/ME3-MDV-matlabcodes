%=============================================================
% Fast Fourier Transform Demo of Leakage:
%  
% Frederic Cegla, 10/01/12
%=============================================================

% Start by defining signals

t=linspace(0,1-1/512,512); % time vector
y=sin(2*pi*10*t);  % 10 Hz signal periodic in time window
y1=sin(2*pi*15.5*t);  % 15.5 Hz signal not periodic in time window

plot(t,y,'r')
hold
plot(t,y1,'b')

legend('10 Hz','15.5Hz')
xlabel('Time (sec)')
ylabel('Amplitude') % plot both


%%
%==================================================================
% Now FFT, and display the different components magnitude and phase
%=================================================================

%take ffts
Fy=fft(y);
Fy1=fft(y1);

T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

figure
% plot magnitude in top plot [ mag(Xk) = sqrt(ak^2+bk^2)]
subplot(2,1,1)
plot(F,abs(Fy1),'-+b')
hold
plot(F,abs(Fy),'-+r')
%plot(F,abs(Fy),'g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Leakage on non-periodic window')
legend('15.5 Hz','10 Hz')
xlim([0 50])
%plot phase in bottom plot  [phase(Xk)= atan(bk/ak) or use Matlab angle function]
subplot(2,1,2)
plot(F,angle(Fy1),'-+b')
hold
plot(F,angle(Fy),'-+r')
%plot(F,angle(Fy),'g')
xlabel('Frequency')
ylabel('Phase (rads)')
legend('15.5 Hz','10 Hz')
xlim([0 50])

%%
%=========================================================================
% Apply hanning Window, FFT, and display the different components magnitude
% and phase
%=========================================================================

%take ffts
Fy=fft(hann(length(y)).*y');
Fy1=fft(hann(length(y1)).*y1');

T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

figure
% plot magnitude in top plot [ mag(Xk) = sqrt(ak^2+bk^2)]
subplot(2,1,1)
plot(F,abs(Fy1),'-+b')
hold
plot(F,abs(Fy),'-+r')
%plot(F,abs(Fy),'g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Minimised leakage by applying a hanning window')
legend('15.5 Hz','10 Hz')
xlim([0 50])
%plot phase in bottom plot  [phase(Xk)= atan(bk/ak) or use Matlab angle function]
subplot(2,1,2)
plot(F,angle(Fy1),'-+b')
hold
plot(F,angle(Fy),'-+r')
%plot(F,angle(Fy),'g')
xlabel('Frequency')
ylabel('Phase (rads)')
legend('15.5 Hz','10 Hz')
xlim([0 50])
%%
