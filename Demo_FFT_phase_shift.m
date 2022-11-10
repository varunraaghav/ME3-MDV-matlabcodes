% Fast Fourier Transform phase shift demo:
% Taking fft of a sine wave and display magnitude and phase
% then change phase and show the effect
% Frederic Cegla, 14/09/19
%=============================================================
% Start by defining signals
t=linspace(0,1-1/512,512); % time vector
y=cos(2*pi*20*t);  % 20 Hz signal
y1=cos(2*pi*20*t-pi/8);  % 20 Hz signal
y2=cos(2*pi*20*t-pi/4);  % 20 Hz signal
y3=cos(2*pi*20*t-3*pi/8);  % 20 Hz signal
y4=cos(2*pi*20*t-pi/2);  % 20 Hz signal

figure
plot(t,y,'b')
hold
plot(t,y1,'r')
plot(t,y2,'g')
plot(t,y3,'m')
plot(t,y4,'y')

legend('20 Hz, \phi=0','20 Hz, \phi=\pi/8', '20 Hz, \phi=\pi/4', '20 Hz, \phi=3\pi/8', '20 Hz, \phi=\pi/2')
xlabel('Time (sec)')
ylabel('Amplitude') % plot 20 and 30 Hz individuall in top plot


%%
% Now FFT, and display the different components magnitude and phase
%%

%take ffts
Fy=fft(y);
Fy1=fft(y1);
Fy2=fft(y2);
Fy3=fft(y3);
Fy4=fft(y4);
T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

figure
% plot magnitude in top plot [ mag(Xk) = sqrt(ak^2+bk^2)]
subplot(2,1,1)
plot(F,abs(Fy),'-+b')
hold
plot(F,abs(Fy1),'-+r')
plot(F,abs(Fy2),'-+g')
plot(F,abs(Fy3),'-+m')
plot(F,abs(Fy4),'-+y')

xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 256]) % display up to Nyquist Frequency
%legend('20 Hz','30 Hz','20 + 30 Hz')
legend('20 Hz, \phi=0','20 Hz, \phi=\pi/8', '20 Hz, \phi=\pi/4', '20 Hz, \phi=3\pi/4', '20 Hz, \phi=\pi/2')

%plot phase in bottom plot  [phase(Xk)= atan(bk/ak) or use Matlab angle function]
subplot(2,1,2)
plot(F,angle(Fy),'-+b')
hold
plot(F,angle(Fy1),'-+r')
plot(F,angle(Fy2),'-+g')
plot(F,angle(Fy3),'-+m')
plot(F,angle(Fy4),'-+y')

%extract phase at 20 Hz
phase_res=[angle(Fy(21)) angle(Fy1(21)) angle(Fy2(21)) angle(Fy3(21)) angle(Fy4(21))]; 

% label plot
xlabel('Frequency')
ylabel('Phase (rads)')
xlim([0 256]) % display up to Nyquist Frequency
legend('20 Hz, \phi=0','20 Hz, \phi=\pi/8', '20 Hz, \phi=\pi/4', '20 Hz, \phi=3\pi/8', '20 Hz, \phi=\pi/2')

% plot phase at 20Hz only 
figure
plot(phase_res)
hold
plot(-[0 pi/8 pi/4 3*pi/8 pi/2],'+r')
xlabel('curve #');
ylabel('Phase (rad)');

% plot a_k and b_k
figure
a_k=[real(Fy(21)) real(Fy1(21)) real(Fy2(21)) real(Fy3(21)) real(Fy4(21))]; 
b_k=[imag(Fy(21)) imag(Fy1(21)) imag(Fy2(21)) imag(Fy3(21)) imag(Fy4(21))]; 
plot([0 pi/8 pi/4 3*pi/8 pi/2], a_k,'+b');
hold
plot([0 pi/8 pi/4 3*pi/8 pi/2], b_k,'+r');
xlabel('Phase (rad)');
ylabel('a_k and b_k');
legend('a_k','b_k')