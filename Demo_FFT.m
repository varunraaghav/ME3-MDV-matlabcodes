% Fast Fourier Transform Application demo 1:
% Taking fft of 2 mixed sine waves and display magnitude and phase 
% Frederic Cegla, 10/01/13
%=============================================================
% Start by defining signals
t=linspace(0,1-1/512,512); % time vector
y1=sin(2*pi*20*t);  % 20 Hz signal
y2=sin(2*pi*30*t);  % 30 Hz signal
y=y1+y2;            % 20 + 30 Hz signal

subplot(2,1,1);     % split plot top plot
plot(t,y1,'b')
hold
plot(t,y2,'r')
legend('20 Hz','30 Hz')
xlabel('Time (sec)')
ylabel('Amplitude') % plot 20 and 30 Hz individuall in top plot

subplot(2,1,2); % bottom plot
plot(t,y,'-b','Linewidth',3) % plot joint 20 + 30 Hz signal in bottom plot
axis([0 1 -2 2])
xlabel('Time (sec)')
ylabel('Amplitude')
legend('20 Hz + 30 Hz')

%%
% Now FFT, and display the different components magnitude and phase
%%

%take ffts
Fy1=fft(y1);
Fy2=fft(y2);
Fy=fft(y);

T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

figure
% plot magnitude in top plot [ mag(Xk) = sqrt(ak^2+bk^2)]
subplot(2,1,1)
plot(F,abs(Fy1),'-+b')
hold
%plot(F,sqrt(real(Fy1).^2+imag(Fy1).^2),'g')
plot(F,abs(Fy2),'-+r')
plot(F,abs(Fy),'-+g')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 256]) % display up to Nyquist Frequency
legend('20 Hz','30 Hz','20 + 30 Hz')
%plot phase in bottom plot  [phase(Xk)= atan(bk/ak) or use Matlab angle function]
subplot(2,1,2)
plot(F,angle(Fy1),'-+b')
hold
%plot(F,atan2(imag(Fy1),real(Fy1)),'m')
plot(F,angle(Fy2),'-+r')
plot(F,angle(Fy),'-+g')
xlabel('Frequency')
ylabel('Phase (rads)')
xlim([0 256]) % display up to Nyquist Frequency
legend('20 Hz','30 Hz','20 + 30 Hz')
%%
% Now reconstruct only the 20 Hz component
%%

res=exp(i*2*pi*20*t).*Fy1(21); %reconstruct time signal from Fy1 component at 20 Hz
figure
plot(t,res/256,'r') %plot reconstructed signal
hold
plot(t,y1,'b--') % plot orignial signal

legend('reconstructed 20Hz signal', 'original 20 Hz signal')
