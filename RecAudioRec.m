% function that records an audio signal when called by button
% press from AudioRecord_and_SpectrumDemo.m
% F. Cegla, 10/01/12
%=========================================================================

function RecAudioRec(source,eventdata)
global y
%sound(y,44100);
%disp('Hello')
t = audiorecorder(44100, 16, 1);
record(t);     % speak into microphone...
pause(2);
stop(t);
y = getaudiodata(t);


%% plot the results
% create split plot
subplot(2,1,1)

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

end
