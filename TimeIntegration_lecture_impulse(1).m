clear all
close all
clc

% System properties
m = 1 ;
k = 1 ;
c = 0.05 ;
w2 = k/m ;
zeta = 0.5*c/sqrt(w2) ;

% State matrix
A = [0, 1; -k/m, -c/m] ;

% Excitation ampliude
Amp = 1 ;
B = [0; 1/m] ;

% Initial conditions
u0 = [0; 0] ;

%% Time-domain parameters
h = 1e-2 ;
tend = 10000 ;
nsamples = 2^nextpow2((tend-h)/h) ;
time = h*linspace(0, nsamples-1, nsamples) ;

fs = 1/h ;
freq = 0.5*fs*linspace(0,1,nsamples/2+1) ;
w = 2*pi*freq ; 

%% The excitation
T0 = 1 ;
W0 = 2*pi/T0 ;

% Time domain excitation
F = Amp*sin(W0*time) ;
nsamples_half_period = 0.5*T0/h ;
F(nsamples_half_period+1:end) = 0 ;
figure ; plot(time, F, 'k', 'linewidth' ,2.0) ;
xlabel('Time (s)') ; ylabel('F(t)') ; xlim([0 2]) ;

% Frequency domain 
Fw = fft(F)/nsamples ;
FwAbs = abs(Fw(1:nsamples/2+1)) ;
FwAbs(2:end-1) = 2*FwAbs(2:end-1) ;
figure ; semilogy(w, FwAbs, 'k', 'linewidth', 2) ;
xlim([0 50]) ;
xlabel('Frequency (rad/s)') ; ylabel('Amplitude |F(\omega)|') ;

Ft_rec = ifft(Fw)*nsamples ;
figure ; plot(time, F, time, Ft_rec, '--') ;
xlim([0 3]) ;
%% Response to half-cosine excitation - Frequency domain
wfft = [w, w(end-1:-1:2)] ;
Xw = Fw./(k+1i*wfft.*c-wfft.^2*m) ;
figure ; semilogy(wfft(1:nsamples/2+1), abs(Xw(1:nsamples/2+1))) ;
xlim([0 50]) ;

Xt_rec = ifft(Xw,'symmetric')*nsamples ;
figure ; plot(time, Xt_rec) ;
xlim([0 50]) ;

%% Response to Half cosine excitation - Time domain
u = zeros(2, length(time)) ;
u(:,1) = u0 ;
% Time integration
for it=1:length(time)-1
    u(:,it+1) = (eye(2)-0.5*h*A)\((eye(2)+0.5*h*A)*u(:,it)+0.5*h*B*(F(it)+F(it+1))) ;
end
figure ; plot(time, u(1,:), 'k', 'linewidth' ,2.0) ;
hold on ; plot(time, Xt_rec, '--', 'linewidth' ,2.0) ;
xlabel('Time (s)') ; ylabel('Displacement') ; %xlim([0 2]) ;
xlim([0 100]) ;

Xwt = fft(u(1,:))/nsamples ;
XwAbs = abs(Xwt(1:nsamples/2+1)) ;
XwAbs(2:end-1) = 2*XwAbs(2:end-1) ;
figure ; semilogy(wfft, abs(Xwt), wfft, abs(Xw), '--') ;
xlim([0 50]) ;

%% Response to half-cosine excitation - Frequency domain Analytical
Fwa = (-(W0*(exp(-(pi*wfft*1i)/W0) + 1))./((wfft + W0).*(wfft - W0))) ;
Fta_rec = (1/h)*ifft(Fwa, 'symmetric');
figure ; plot(time, F, time, Fta_rec, '--') ;
xlim([0 3]) ;

Xwta = Fwa./(k+1i*wfft.*c-wfft.^2*m)/nsamples/h ;
figure ; semilogy(wfft, abs(Xwta), 'k', wfft, abs(Xwt), '--', 'linewidth' ,2.0) ;
xlim([0 50]) ;
xlabel('Frequency (rad/s)') ; ylabel('Displacement amplitude |X(\omega)|') ;

Xta_rec = ifft(Xwta,'symmetric')*nsamples ;
figure ; plot(time, u(1,:), 'k', time(1:40:end), Xta_rec(1:40:end), '*', 'linewidth' ,2.0) ;
xlim([0 100]) ;
xlabel('Time (s)') ; ylabel('Displacement') ; 
