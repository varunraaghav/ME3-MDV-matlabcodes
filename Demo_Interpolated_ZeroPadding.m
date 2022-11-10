% Demo that illustrates zero padding of the frequency domain vector as 
% tool to interpolate time domain data
% Frederic Cegla, 10/09/19
%=============================================================

% Start by defining signals
t=linspace(0,1-1/512,512); % time vector
%t=[0:1/512:1]
y=sin(2*pi*20*t).*hann(length(t))';

%%
% Now FFT, and display the different components magnitude and phase
%%

%take ffts
Fy=fft(y);
T=t(end); % period
F=[0:1/T:(length(t)-1)/T]; %define frequency vector

figure
subplot(2,1,1)
plot(t,y,'b','Linewidth',3)
%hold
%plot(t1,y1,'r')
%legend('1 sec.')
xlabel('Time (sec)')
ylabel('Amplitude') % plot signals
subplot(2,1,2)
plot(F(1:end/2),abs(Fy(1:end/2)),'b','Linewidth',3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% produce zero padded frequency vector
N_original=length(Fy);
Fy2=zeros(1,8*N_original);
Fy2(2:N_original/2)=8.*Fy(2:N_original/2);
Fy2(end-(N_original/2):end)=8.*Fy(end-(N_original/2):end);
df=1;
F2=[0:df:(8*N_original-1)*df];
y2=real(ifft(Fy2));
Fmax=F2(end)+df;
t2=[0:1/Fmax:(length(F2)-1)/Fmax];

figure
plot(t,y,'-b+')
hold
plot(t2,y2,'-r')

%% do it for a picture

%read in picture
pic=imread('Horse_detail.jpg');
%convert to grayscale
I = rgb2gray(pic(1:256,1:256,:));

%plot original picture
figure
imagesc(I)

%compute fft
F_im=fft2(I);

%retrieve size of fft
N_im=length(F_im(:,1));

%create new empty increased fft vector
upscale_factor=8;
F_new=zeros(N_im*upscale_factor,N_im*upscale_factor);
%populate the four quadrants with the existing data
F_new(1:N_im/2,1:N_im/2)=upscale_factor^2*F_im(1:N_im/2,1:N_im/2);
F_new(end-N_im/2:end,end-N_im/2:end)=upscale_factor^2*F_im(end-N_im/2:end,end-N_im/2:end);
F_new(1:N_im/2,end-N_im/2:end)=upscale_factor^2*F_im(1:N_im/2,end-N_im/2:end);
F_new(end-N_im/2:end,1:N_im/2)=upscale_factor^2*F_im(end-N_im/2:end,1:N_im/2);
%calculated inverse fft
New_I=floor(real(ifft2(F_new)));
% display upsampled picture
figure
imagesc(New_I)
title('with interpolation')

%without interpol
IBig=zeros(size(I).*upscale_factor);
dim_pic=size(I);
for count=1:dim_pic(1);
    for count1=1:dim_pic(2);
        for count2=1:8;
            for count3=1:8
                IBig(((count-1)*upscale_factor+count2),((count1-1)*upscale_factor)+count3)=I(count,count1);     
            end
        end
    end
end
figure
imagesc(IBig)
title('without interpolation')


