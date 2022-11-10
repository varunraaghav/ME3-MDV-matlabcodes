clear all
close all
clc

%% System's definition
M = eye(2) ;
K = [2, -1; -1, 2] ;
[X, w2] = eig(K, M, 'vector') ;

%% Orthogonality properties
X'*M*X
X'*K*X

%% Proportional Damping
a = 0.01 ;
b = 0 ;
C = a*K+b*M ;
zeta = 0.5*diag(X'*C*X)./sqrt(w2)

%% FRFs - Physical
w = linspace(0, 2*sqrt(w2(2)), 1e4) ;
Hw = zeros(2,2,length(w)) ;
Hw_modal = zeros(2,2,length(w)) ;
Hmag = zeros(2,2,length(w)) ;
for iw = 1:length(w)
    for iforce = 1:2
        F = zeros(2,1) ;
        F(iforce) = 1 ;
        Hw(:,iforce,iw) = (K+1i.*w(iw)*C-w(iw)^2*M)\ F ;
        Hw_modal(:,iforce,iw) = (X'*(K+1i.*w(iw)*C-w(iw)^2*M)*X)\(X'*F) ;
    end
end

% Display
figure ; semilogy(w, squeeze(abs(Hw(1,1,:))), 'k', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Amplitude') ;
box on ; xlim([w(1) w(end)]) ;

figure ; plot(w, squeeze(angle(Hw(1,1,:))), 'k', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Phase (rad)') ;
box on ; xlim([w(1) w(end)]) ;

figure ; plot(squeeze(real(Hw(1,1,:))), squeeze(imag(Hw(1,1,:))), 'k', 'linewidth', 2) ;
xlabel('Real part') ; ylabel('Imag part') ; box on ;

figure ; semilogy(w, squeeze(abs(Hw(2,1,:))), 'k', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Amplitude') ;
box on ; xlim([w(1) w(end)]) ;

figure ; plot(w, (squeeze(angle(Hw(2,1,:)))), 'k', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Phase (rad)') ;
box on ; xlim([w(1) w(end)]) ;

figure ; plot(squeeze(real(Hw(2,1,:))), squeeze(imag(Hw(2,1,:))), 'k', 'linewidth', 2) ;
xlabel('Real part') ; ylabel('Imag part') ; box on ;

x_f1_1 = (X(1,:)'*ones(1,length(w))).*squeeze(Hw_modal(:,1,:)) ;
figure ; semilogy(w, squeeze(abs(Hw(1,1,:))), 'k', 'linewidth', 2) ;
hold on ; semilogy(w, abs(x_f1_1(1,:)), '--', 'linewidth', 2) ;
hold on ; semilogy(w, abs(x_f1_1(2,:)), '--', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Amplitude') ;
box on ; xlim([w(1) w(end)]) ;

figure ; plot(w, squeeze(real(Hw(1,1,:))), 'k', 'linewidth', 2) ;
hold on ; semilogy(w, real(x_f1_1(1,:)), '--', 'linewidth', 2) ;
hold on ; semilogy(w, real(x_f1_1(2,:)), '--', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Real part') ;
hold on ; plot([w(1) w(end)], [0 0], '--k') ;
box on ; xlim([0.5 2]) ;

figure ; plot(w, squeeze(imag(Hw(1,1,:))), 'k', 'linewidth', 2) ;
hold on ; semilogy(w, imag(x_f1_1(1,:)), '--', 'linewidth', 2) ;
hold on ; semilogy(w, imag(x_f1_1(2,:)), '--', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Imag part') ;
hold on ; plot([w(1) w(end)], [0 0], '--k') ;
box on ; xlim([0.5 2]) ;

x_f1_2 = (X(2,:)'*ones(1,length(w))).*squeeze(Hw_modal(:,1,:)) ;
figure ; semilogy(w, squeeze(abs(Hw(2,1,:))), 'k', 'linewidth', 2) ;
hold on ; semilogy(w, abs(x_f1_2(1,:)), '--', 'linewidth', 2) ;
hold on ; semilogy(w, abs(x_f1_2(2,:)), '--', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Amplitude') ;
box on ; xlim([w(1) w(end)]) ;

figure ; plot(w, squeeze(real(Hw(2,1,:))), 'k', 'linewidth', 2) ;
hold on ; semilogy(w, real(x_f1_2(1,:)), '--', 'linewidth', 2) ;
hold on ; semilogy(w, real(x_f1_2(2,:)), '--', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Real part') ;
hold on ; plot([w(1) w(end)], [0 0], '--k') ;
box on ; xlim([0.5 2]) ;

figure ; subplot(1,2,1) ; semilogy(w, squeeze(abs(Hw(2,1,:))), 'k', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Amplitude') ;
box on ; xlim([w(1) w(end)]) ;
subplot(1,2,2) ; semilogy(w, squeeze(abs(Hw(1,2,:))), 'k', 'linewidth', 2) ;
xlabel('Frequency') ; ylabel('Amplitude') ;
box on ; xlim([w(1) w(end)]) ;

