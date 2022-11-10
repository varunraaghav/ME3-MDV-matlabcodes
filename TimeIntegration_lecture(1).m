clear all
close all
clc

%%
w = 2*pi ;
A = [0, 1; -w^2, 0] ;

u0 = [1; 0] ;

%% Reference
[t_ref, u_ref] = ode45(@(t,x) A*x, [0 5], u0) ;

%% Euler Forward
h = [5e-3, 1e-2] ;

figure ;
for ih = 1:length(h)
    time = 0:h(ih):5 ;
    u = zeros(2, length(time)) ;
    u(:,1) = u0 ;
    for it=1:length(time)-1
        u(:,it+1) = u(:,it)+h(ih)*A*u(:,it) ;
    end
    plot(time, u(1,:), 'linewidth', 2.0) ; hold on ;
end
plot(time, cos(w*time), 'k', 'linewidth', 2.0) ; hold on ;
xlabel('Time (s)') ; ylabel('Displacement') ;
% plot([1 1], [min(u(1,:)) max(u(1,:))], 'k') ; 

figure ;
plot(cos(linspace(0, 1.999*pi, 1e3)), sin(linspace(0, 1.999*pi, 1e3)),'k', 'linewidth', 2) ; hold on
for ih = 1:length(h)
    [~,lambda] = eig(eye(2)+h(ih)*A,'vector') ;
    plot(real(lambda), imag(lambda), '*', 'linewidth', 2.0, 'Markersize', 12) ;
end

theta = linspace(0, 1.999*pi, 1e3) ;
figure ; plot([0 0], [0 0], '--k') ; hold on ;
xlim([-3 3]) ; ylim([-3 3]) ; axis equal ; box on ;
fill(cos(theta)-1, sin(theta), [0.85,0.33,0.10], 'edgecolor', 'none', 'facealpha', 0.5) ; hold on ;
plot([0 0], [-4 4], '--k') ; plot([-4 4], [0 0], '--k') ; 
plot(cos(theta)-1, sin(theta), 'k', 'linewidth', 2.) ;
xlabel('Real(h \lambda)') ; ylabel('Im(h \lambda)') ;

%% Euler Backward
h = [5e-3, 1e-2] ;

figure ;
for ih = 1:length(h)
    time = 0:h(ih):5 ;
    u = zeros(2, length(time)) ;
    u(:,1) = u0 ;
    for it=1:length(time)-1
        u(:,it+1) = (eye(2)-h(ih)*A)\u(:,it) ;
    end
    plot(time, u(1,:), 'linewidth', 2.0) ; hold on ;
end
plot(time, cos(w*time), 'k', 'linewidth', 2.0) ; hold on ;
xlabel('Time (s)') ; ylabel('Displacement') ;

figure ;
plot(cos(linspace(0, 1.999*pi, 1e3)), sin(linspace(0, 1.999*pi, 1e3)),'k', 'linewidth', 2) ; hold on
for ih = 1:length(h)
    [~,lambda] = eig((eye(2)-h(ih)*A)\eye(2),'vector') ;
    plot(real(lambda), imag(lambda), '*', 'linewidth', 2.0, 'Markersize', 12) ;
end

theta = linspace(0, 1.999*pi, 1e3) ;
figure ; plot([0 0], [0 0], '--k') ; hold on ;
xlim([-3 3]) ; ylim([-3 3]) ; axis equal ; box on ;
xrange = get(gca, 'Xlim') ; 
fill([xrange(1) xrange(1) xrange(2) xrange(2)], [-3 3 3 -3], [0.85,0.33,0.10], 'edgecolor', 'none', 'facealpha', 0.5) ; hold on ;
plot([0 0], [-4 4], '--k') ; plot([-4 4], [0 0], '--k') ; 
plot(1-cos(theta), -sin(theta), 'k', 'linewidth', 2.) ;
fill(1-cos(theta), -sin(theta), [1,1,1], 'edgecolor', 'none') ;
xlabel('Real(h \lambda)') ; ylabel('Im(h \lambda)') ;

%% Trapezoidal rule
h = [5e-3, 1e-2] ;

figure ; 
plot(time, cos(w*time), 'k', 'linewidth', 2.0) ; hold on ;
for ih = 1:length(h)
    time = 0:h(ih):5 ;
    u = zeros(2, length(time)) ;
    u(:,1) = u0 ;
    for it=1:length(time)-1
        u(:,it+1) = (eye(2)-0.5*h(ih)*A)\(eye(2)+0.5*h(ih)*A)*u(:,it) ;
    end
    plot(time, u(1,:), '--', 'linewidth', 2.0) ; hold on ;
end
xlabel('Time (s)') ; ylabel('Displacement') ;
% plot([1 1], [min(u(1,:)) max(u(1,:))], 'k') ; 

figure ;
plot(cos(linspace(0, 1.999*pi, 1e3)), sin(linspace(0, 1.999*pi, 1e3)),'k', 'linewidth', 2) ; hold on
for ih = 1:length(h)
    [~,lambda] = eig((eye(2)-0.5*h(ih)*A)\(eye(2)+0.5*h(ih)*A),'vector') ;
    plot(real(lambda), imag(lambda), '*', 'linewidth', 2.0, 'Markersize', 12) ;
end

figure ; plot([0 0], [0 0], '--k') ; hold on ;
xlim([-3 3]) ; ylim([-3 3]) ; axis equal ; box on ;
xrange = get(gca, 'Xlim') ; 
fill([xrange(1) xrange(1) 0 0], [-3 3 3 -3], [0.85,0.33,0.10], 'edgecolor', 'none', 'facealpha', 0.5) ; hold on ;
plot([0 0], [-4 4], '--k') ; plot([-4 4], [0 0], '--k') ; 
xlabel('Real(h \lambda)') ; ylabel('Im(h \lambda)') ;
