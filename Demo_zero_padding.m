%=============================================================
% Fast Fourier Transform Demo of zero padding:
%  
% Frederic Cegla, 23/08/19
%=============================================================

t=linspace(0,0.05-0.05/512,512); % time vector
% maximum frequency in spectrum Fmax=1/dt=1/(0.05/512)=10240 Hz
% however maximum frequency that can be resolved
% is the Nyquist frequency Fn=Fmax/2=10240/2=5120 Hz
% frequency resolution df=1/T=1/0.05=20Hz

T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

% Define one signal with frequency 550 Hz
y=sin(2*pi*550*t);  % 50 Hz signal periodic in time window

%Define one signal with frequency 560 Hz
y1=sin(2*pi*560*t);  % 400 Hz signal periodic in time window

% plot and check signals if desired 
check=0;
if check==1
    plot(t,y,'b')
    hold
    plot(t,y1,'r')
    legend('550 Hz','560Hz')
    xlabel('Time (sec)')
    ylabel('Amplitude') 
end

%combine both signals
y=y+y1;

plot(t,y,'b','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('Time (sec)','Fontsize',16)
ylabel('Amplitude','Fontsize',16)
title('550 and 560 Hz mixed sine waves','Fontsize',16)

%%
%==================================================================
% Now FFT, and display the different components magnitude and phase
%=================================================================

%take ffts
Fy=fft(y);
figure
plot(F,abs(Fy),'b-+','Linewidth',2);
set(gca,'Fontsize',16)
xlabel('Frequency (Hz)','Fontsize',16)
ylabel('Amplitude','Fontsize',16)
xlim([0 1000])

% sample for longer

t=linspace(0,0.4-0.4/4096,4096); % time vector
% maximum frequency in spectrum Fmax=1/dt=1/(0.4/4096)=10240 Hz
% however maximum frequency that can be resolved
% is the Nyquist frequency Fn=Fmax/2=10240/2=5120 Hz
% frequency resolution df=1/T=1/0.4=2.5Hz

T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

% Define one signal with frequency 550 Hz
y=sin(2*pi*550*t);  % 50 Hz signal periodic in time window

%Define one signal with frequency 560 Hz
y1=sin(2*pi*560*t);  % 400 Hz signal periodic in time window
y=y+y1;
plot(t,y,'b','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('Time (sec)','Fontsize',16)
ylabel('Amplitude','Fontsize',16)
title('550 and 560 Hz mixed sine waves','Fontsize',16)
%take ffts
Fy=fft(y);
figure
plot(F,abs(Fy),'b-+','Linewidth',2);
set(gca,'Fontsize',16)
xlabel('Frequency (Hz)','Fontsize',16)
ylabel('Amplitude','Fontsize',16)
xlim([380 720])

% zero padded result
y=zeros(1,length(y));
y(1:512)=sin(2.*pi.*550.*t(1:512))+sin(2.*pi.*560.*t(1:512));
plot(t,y,'b','Linewidth',2)
set(gca,'Fontsize',16)
xlabel('Time (sec)','Fontsize',16)
ylabel('Amplitude','Fontsize',16)
title('550 and 560 Hz mixed sine waves zero padded','Fontsize',16)
%take ffts
Fy=fft(y);
figure
plot(F,abs(Fy),'b-+','Linewidth',2);
set(gca,'Fontsize',16)
xlabel('Frequency (Hz)','Fontsize',16)
ylabel('Amplitude','Fontsize',16)
xlim([380 720])
