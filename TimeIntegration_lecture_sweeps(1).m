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

% Road periodic profile
Amp = 1 ;
B = [0; 1/m] ;
w = linspace(0, 5, 1e3) ;

%% Frequency domain solution
Xfreq = Amp*1./(k-w.^2*m+1i*w.*c) ;

%% Time-domain solution
u0 = [abs(Xfreq(1)); 0] ;
h = 1e-2 ;
Xtime = zeros(1, length(w)) ;
for iw = 2:length(w)
    T = 2*pi/w(iw) ;
    time = 0:h:5*T ; % simulate 5 oscillation periods
    nsamples_period = ceil(T/h) ;
    u = zeros(2, length(time)) ;
    u(:,1) = u0 ;
    F = Amp*cos(w(iw)*time) ;
    for it=1:length(time)-1
        u(:,it+1) = (eye(2)-0.5*h*A)\((eye(2)+0.5*h*A)*u(:,it)+0.5*h*B*(F(it)+F(it+1))) ;
    end
    fode = @(t, z) [A*z+B*Amp*cos(w(iw)*t)] ;
    [~, xref] = ode45(fode, time, u0) ;
    Xtime(1,iw) = max(abs(u(1,end-nsamples_period+1:end))) ;
    u0 = u(:,end) ;
end

figure ; semilogy(w, abs(Xfreq), 'k', 'linewidth', 2) ;
hold on ; semilogy(w, Xtime, '--', 'linewidth', 2) ;
xlabel('\omega') ; ylabel('Response amplitude') ;

%% Response to sine sweeps with different rates
u0 = [0; 0] ;
h = 1e-2 ;
rates = 0.1*[0.1 1 3 5]/2/pi ;
f0 = 0.1/2/pi ;
fend = 5/2/pi ;
figure ;
for ir = 1:length(rates)
    tend =  (fend-f0)*(60/rates(ir)) ;
    k = sign(rates(ir))*abs(fend-f0)/tend ;
    time = 0:h:tend ;
    finst = k*time + f0 ;
    psi = 2*pi*(f0.*time+0.5*k*time.^2) ;
    
    u = zeros(2, length(time)) ;
    u(:,1) = u0 ;
    F = Amp*sin(psi) ;
    for it=1:length(time)-1
        u(:,it+1) = (eye(2)-0.5*h*A)\((eye(2)+0.5*h*A)*u(:,it)+0.5*h*B*(F(it)+F(it+1))) ;
    end
    plot(2*pi*finst, u(1,:), 'linewidth', 2) ; hold on ; drawnow ;
end
hold on ; plot(w, abs(Xfreq),  'k', 'linewidth', 2) ;
xlim([2*pi*f0 3]) ; xlabel('Instantaneous frequency (rad/s)') ; ylabel('Displacement') ;

figure ; plot(time, F, 'k', 'linewidth' ,2.0) ;
xlabel('Time (s)') ; ylabel('F(t)') ; xlim([0 150]) ;

