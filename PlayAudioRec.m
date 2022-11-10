% function that plays the recorded audio object r when called by button
% press from AudioRecord_and_SpectrumDemo.m
% F. Cegla, 10/01/12
%=========================================================================

function PlayAudioRec(source,eventdata)
global y
sound(y,44100);
disp('Hello')
end
