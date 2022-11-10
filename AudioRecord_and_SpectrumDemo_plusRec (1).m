% =======================================================================
% Demo Program that records an audio signal and displays its spectrum
% F.Cegla, Department of Mechanical Engineering
% Imperial College London
% January 2012
% =======================================================================
% for the sound to be playable by another function the variabl y that contains 
% information needs to
% be made global so that the PlayAudioRec function can access it
global y
figure
%create audio object
r = audiorecorder(44100, 16, 1);
% create split plot
subplot(2,1,1)
record(r);     % speak into microphone...
pause(2);
stop(r);
% get recorded audio data and plot
y = getaudiodata(r);
t=[0:1/44100:(length(y)-1)/44100];
plot(t,y)
xlabel('Time (Seconds)'); ylabel('Amplitude (V)');
title('Recorded time domain signal')

%do FFT analysis and display spectrum in second plot
subplot(2,1,2)
Fy=fft(hann(length(y)).*y);
T=t(end);
%create frequency vector from time information
F=[0:1/T:(length(y)-1)/T];
% plot spectrum up to the nyquist frequency = Fmax/2
plot(F(1:floor(length(F)/2)),abs(Fy(1:floor(length(F)/2))))
xlabel('Frequency(Hz)'); ylabel('Amplitude');
title('Spectrum of recorded signal')
xlim([0 5000])
zoom on
%%
% include button to re-play recorded signal

% then define button, which when pressed will call PlayAudioRec and play
% back the recorded signal
pbh1 = uicontrol('Style','pushbutton','String','Play',...
                'Units','normalized',...
                'Position',[.92 .7 .05 .1],'Callback',@PlayAudioRec);
            
%%
% include button to record signal

% then define button, which when pressed will call PlayAudioRec and play
% back the recorded signal
pbh1 = uicontrol('Style','pushbutton','String','Record',...
                'Units','normalized',...
                'Position',[.92 .3 .05 .1],'Callback',@RecAudioRec);