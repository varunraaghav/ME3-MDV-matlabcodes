% Demo that illustrates response of an SDOF to a user defined impulse
% MDV course demo
% Frederic Cegla, 16/03/21
%=============================================================

% Start by defining signals
t=linspace(0,1-1/512,512); % time vector
y=zeros(size(t));
p=t(1:5);
y(51:55)=ones(1,length(p)).*hann(length(p))';
% p=t(1:25);
% y(51:75)=ones(1,length(p)).*hann(length(p))';

%y(101:105)=-ones(1,length(p)).*hann(length(p))';
y1=y;

%%
% Now FFT, and display the different components magnitude and phase
%%
%take fft of input signal
Fy1=fft(y1);
% compute frequncy consitutents
T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector
%% plot
figure
subplot(3,1,1)
plot(t,y1,'b','Linewidth',3)
%legend('1 sec.')
xlabel('Time (sec)')
ylabel('Amplitude') % plot signals
subplot(3,1,2)
plot(F(1:end/2),abs(Fy1(1:end/2)),'b','Linewidth',3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% Response FRF for an SDOF
m1=0.05; % mass    in kg
f1= 100;  %hz
c1=0.1 ;  % viscous damping ratio 

w=2*pi*F; % frequency vector (rad /sec)
w1=2*pi*f1; % resonant frequency, rad/sec
k1=m1*w1*w1; % hence stiffness k1 in N/m
x=w/w1;    %freq ratio
a=(1-x.*x); % real denominator
b=2.*c1.*x ;  % imag denominator
denom = complex(a,b);
% compute FRF   
h=1/k1./denom;
subplot(3,1,3)
plot(F(1:end/2),abs(h(1:end/2)),'b','Linewidth',3)
xlabel('Frequency (Hz)')
ylabel('SDOF response Amplitude')

%% Multipy FRF with excitation spectrum
FResponse1=h.*Fy1;
% only the spectrum up to the Nyquist frequency is correctly computed
% produce spectrum by using mirror of Fourrier coefficients below the 
% Nyquist frequency  
FResponse1(end/2+1:end)=conj(FResponse1(end/2:-1:1));
% Calculate inverse FFT to compute combined response
TRes1=real(ifft(FResponse1));
%Plot result
figure
subplot(2,1,1)
plot(t,y1)
xlabel('Time (sec)')
ylabel('Excitation Amplitude') % plot excitation signal
subplot(2,1,2)
plot(t,TRes1)
xlabel('Time (sec)')
ylabel('Response Amplitude') % plot Response signal