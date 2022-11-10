% Demo that shows the bandwidth in the frequency spectrum
% of different hann window shaped impulses
% Demo for the MDV course
% Frederic Cegla, 16/03/21
%=============================================================

% Start by defining signals
t=linspace(0,1-1/512,512); % time vector
%t=[0:1/512:1]
y=zeros(size(t));
p=t(1:5);
y(1:5)=ones(1,length(p)).*hann(length(p))';
y1=y;
p=t(1:10);
y(1:10)=ones(1,length(p)).*hann(length(p))';
y2=y;
p=t(1:20);
y(1:20)=ones(1,length(p)).*hann(length(p))';
y3=y;

%%
% Now FFT, and display the different components magnitude and phase
%%

%take ffts
Fy1=fft(y1);
Fy2=fft(y2);
Fy3=fft(y3);

T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

figure
subplot(2,1,1)
plot(t,y1,'b','Linewidth',3)
hold
plot(t,y2,'r','Linewidth',3)
plot(t,y3,'k','Linewidth',3)
%legend('1 sec.')
xlabel('Time (sec)')
ylabel('Amplitude') % plot signals
subplot(2,1,2)
plot(F(1:end/2),abs(Fy1(1:end/2)),'b','Linewidth',3)
hold
plot(F(1:end/2),abs(Fy2(1:end/2)),'r','Linewidth',3)
plot(F(1:end/2),abs(Fy3(1:end/2)),'k','Linewidth',3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
