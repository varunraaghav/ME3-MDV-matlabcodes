% ===============================================================
% ===============================================================
% ===============================================================
% 
%  ME3 - Machine Dynamics and Vibration - SDOF demo
%
%  frf regeneration  with complex function
%  and display of different forms of FRF and Time domain output
%
%  F. Cegla September 2012
%
% ===============================================================
% ===============================================================
% ===============================================================

f = [1:500] ;  % freq range in Hz 1 Hz increments
m1 = 1 ; % mass    in kg
f1 = 100 ;  % hz
c1 = 0.05 ;  % viscous damping ratio 

w=2*pi*f; % rad /sec
w1=2*pi*f1; % res freq, rad/sec
k1=m1*w1*w1; % hence stiffness k1 in N/m
x=w/w1;    % freq ratio
a=(1-x.*x); % real denominator
b=c1.*x ;  % imag denominator
denom = complex(a,b);
   
h=1/k1./denom;
rrecep=real(h) ;
irecep=imag(h); 
mag=abs(h);
phase =atan2(irecep,rrecep); % correct to fourth quadrant
phase=phase.*180/pi;  %convert to degrees

%% Plot results

% define a plot with three subplot regions and make the top on active
SS=get(0,'Screensize');
figure('Position',SS)
h0a=subplot(2,2,1);
semilogy(f,mag,'r','Linewidth',2) %plot receptance of mass1
xlabel('Frequency Hz','Fontsize',14)
ylabel('Log Recpetance (log (x/F)) ','Fontsize',14)
FRFylim=[min(min(abs(mag))) max(max(abs(mag)))];
%ylim(FRFylim)
grid on
%plot sketch of the system  
subplot(2,2,3)
plot(f,phase,'g','Linewidth',2)
grid on
xlabel('Frequency Hz','Fontsize',14)
ylabel('Phase (deg) ','Fontsize',14)


% plot a nyquist plot
h0c=subplot(2,2,2);
plot (rrecep,irecep,'r+-')
xlabel(' Real Part')
ylabel ('imag part')
title(' Nyquist Plot of SDOF')
grid on 
xlabel('Real axis','Fontsize',14)
ylabel('Imaginary axis ','Fontsize',14)
title('Nyquist plot','Fontsize',14)

% Inform user of who to use the demo
instructions=helpdlg(strcat('Click at a point in the top left graph to select a frequency',...
    ' at which the mass motion is simulated in the bottom plot.',...
    'Wait for the simulation to finish before selecting a new frequency.',...
    'Press any key to exit')); 
waitfor(instructions)

%Radiobutton for selection of acceleration or Displacement option
%==========================================================================
% Create one radio buttons 
% u0 = uicontrol('Style','Radio','String','Acceleration/Force',...
%     'pos',[200 500 150 30],'HandleVisibility','off');

%set initial as displacement FRF
%set(u0,'Value',0)
%==========================================================================
%Loop for interaction with the user
%==========================================================================
w=0; % define entry condition to the while loop
fselection=0; %turn to 1 if frequency has been selected
while w==0
    % catch the location of where the mouse is clicked
    % if a key is pressed exit the program
    w=waitforbuttonpress;
    if w == 0
        disp('Button click')
    else
        disp('Key press')
        close all
        break
    end
    %retrieve button control
%     b1=get(u0,'Value')  % value button 1
%     
%     if b1==1
%             disp('displacement off')
%             set(u0,'Value',1)
%             disp('acceleration on')
%             AccelFRF
%     elseif b1==0
%             disp('displacement on')
%             set(u0,'Value',0)
%             disp('acceleration off')
%     end
    
    
      
    % extract the selected frequency where the mouse was clicked
    frequest=get(h0a,'CurrentPoint');
    fchoice=frequest(1,1);
    if fchoice<f(1)
      fchoice=f(10)
     end
     if fchoice>f(end)
         fchoice=f(end-10)
     end
    frequest(1,1)=fchoice;
     
    omega=2*pi*frequest(1,1);
       
    % indicate selected frequency on FRF and Nyquist plots
    if fselection ==1
        %update lines if they already exist
        set(hl1,'XData',[frequest(1,1) frequest(1,1)]);
        set(hl2,'XData',[0 real(h(floor(frequest(1,1))))]);
        set(hl2,'YData',[0 imag(h(floor(frequest(1,1))))]);
    else
        %create lines in case they do not already exist
        h0a=subplot(2,2,1);
        hl1=line([frequest(1,1) frequest(1,1)],FRFylim,'Color',[0 0 0 ],'Linewidth',2);
        h0c=subplot(2,2,2);
        hl2=line([0 real(h(floor(frequest(1,1))))],[0 imag(h(floor(frequest(1,1))))],'Color',[0 0 0 ],'Linewidth',2);
        fselection=1;
    end
    
    % normalise x
    x=x./max(abs(x));
       
    % set bottom plot to current
    h0b=subplot(2,2,4);
    % up date the selected frequency in the title
    title(strcat('Time signals at:  ',num2str(floor(frequest(1,1))),' Hz'),'Fontsize',16)
    % Plot time histories
    
    t=0:1/frequest(1,1)/20:3/frequest(1,1);
    F=cos(omega.*t); %define Force
    x1=mag(floor(frequest(1,1))).*cos(omega.*t+phase(floor(frequest(1,1)))/180*pi)*1000;
    plot(t,F,'r')
    hold
    plot(t,x1,'g')
    hold
    xlabel('Time (secs)','Fontsize',14)
    ylabel('Amplitude in (N-red and mm - green)','Fontsize',14)
    
end

