clear all
close all
clc

%% System's definition
m1 = 1 ;
m2 = 0.05*m1 ;
k1 = 1 ;
k2 = k1*m2/m1 ;
c1 = 2e-3;
c2 = 0.0 ;
M = [m1,0;0,m2] ;
C = [c1+c2, -c2; -c2, c2] ;
K = [k1+k2, -k2; -k2, k2] ;
[X, w2] = eig(K, M, 'vector') ;

%% Proportional Damping
zeta = 0.5*diag(X'*C*X)./sqrt(w2) ;

%% FRFs - Physical
w = linspace(0, 2*sqrt(w2(2)), 1e4) ;
Htva = zeros(2,length(w)) ;
Hsdof = zeros(1,length(w)) ;
F = zeros(2,1) ;
F(1) = 1 ;
for iw = 1:length(w)
    Hsdof(:,iw) = (k1+1i.*w(iw)*C(1,1)-w(iw)^2*m1)\ F(1) ;
    Htva(:,iw) = (K+1i.*w(iw)*C-w(iw)^2*M)\F ;
end

figure ; semilogy(w, abs(Hsdof), '--k', 'linewidth', 2) ;
hold on ; semilogy(w, abs(Htva(1,:)), 'k', 'linewidth', 2) ;
xlabel('\omega/\omega_n') ; ylabel('|x_1|') ;
box on ; xlim([w(1) w(end)]) ;

figure ; semilogy(w, abs(Htva(2,:)), 'k', 'linewidth', 2) ;
xlabel('\omega/\omega_n') ; ylabel('|x_2|') ;
box on ; xlim([w(1) w(end)]) ;

%% Mass ratio and natural frequencies
mu = linspace(0.0001, 5, 1e5) ;
wp = (1+0.5*mu+sqrt(mu+0.25*mu.^2)) ;
wm = (1+0.5*mu-sqrt(mu+0.25*mu.^2)) ;
figure ; plot(mu, wp, 'k', 'linewidth', 2.0) ;
hold on ; plot(mu, wm, 'k', 'linewidth', 2.0) ;
xlabel('\mu = m_2/m_1') ; ylabel('\omega^2_{1,2}/\omega^2_n') ;

%% Equal-peak method
k2 = 0.0454 ;
c2 = 0.0128 ;
M = [m1,0;0,m2] ;
C = [c1+c2, -c2; -c2, c2] ;
K = [k1+k2, -k2; -k2, k2] ;
w = linspace(0, 2*sqrt(w2(2)), 1e4) ;
Heqp = zeros(2,length(w)) ;
F = zeros(2,1) ;
F(1) = 1 ;
for iw = 1:length(w)
    Heqp(:,iw) = (K+1i.*w(iw)*C-w(iw)^2*M)\F ;
end

figure ; semilogy(w, abs(Hsdof), '--k', 'linewidth', 2) ;
hold on ; semilogy(w, abs(Htva(1,:)), 'k', 'linewidth', 2) ;
semilogy(w, abs(Heqp(1,:)), 'linewidth', 2) ;
xlabel('\omega/\omega_n') ; ylabel('|x_1|') ;
box on ; xlim([w(1) w(end)]) ;

figure ; semilogy(w, abs(Htva(2,:)), 'k', 'linewidth', 2) ;
semilogy(w, abs(Heqp(2,:)), 'linewidth', 2) ;
xlabel('\omega/\omega_n') ; ylabel('|x_2|') ;
box on ; xlim([w(1) w(end)]) ;


