%=============================================================
% Fast Fourier Transform Demo of Aliasing:
%  
% Frederic Cegla, 10/01/12
%=============================================================

t=linspace(0,1-1/512,512); % time vector
% maximum frequency in spectrum Fmax=1/dt=512 Hz
% however maximum frequency that can be resolved
% is the Nyquist frequency Fn=Fmax/2=512/256 Hz

% Define one signal that is resolved at 50 Hz

y=sin(2*pi*50*t);  % 50 Hz signal periodic in time window

%Define one signal that is not resolved at e.g. 400 Hz
y1=sin(2*pi*400*t);  % 400 Hz signal periodic in time window

plot(t,y,'b')
hold
plot(t,y1,'r')

legend('50 Hz','400Hz')
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
plot(F,abs(Fy1),'b')
hold
plot(F,abs(Fy),'r')
plot([256 256],[0 300],'k--')
legend('400 Hz','50 Hz','Nyquist Frequency')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Aliasing')
%xlim([0 256])

%plot phase in bottom plot  [phase(Xk)= atan(bk/ak) or use Matlab angle function]
subplot(2,1,2)
%xlim([0 256])
plot(F,angle(Fy1),'b')
hold
plot(F,angle(Fy),'r')
xlabel('Frequency')
ylabel('Phase (rads)')

