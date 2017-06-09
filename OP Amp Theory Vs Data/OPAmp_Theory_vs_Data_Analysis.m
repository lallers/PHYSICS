%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Physics Junior Lab, 2016                           /
%     University of New Mexico                           /
%/////////////////////////////////////////////////////////

clear all; clc; close all;
%<<<<<<<<<< Functions Used >>>>>>>>>>>>
F.F1 =@(R,C) 1/(2.*log(3).*R.*C) ; %(R3,C) Oscilation frequency determined by components R3 and C in Relaxation Oscilaltor
F.F2 = @(R1,R2,C) 1.44/((R1 + 2*R2)*C); %Used for 55 Timer Chip Circuit Oscillation Frequency                           
F.DC = @(R1,R2)  (R1 + R2)/(R1 + 2*R2); %Duty Cycle for Timer Chip
F.y =@(R,C) 1/(R*C);
F.TP =@(R1,R2,C,Vc) (R1+R2)*C*log((5-.5*Vc)/(5-Vc));
F.T = @(Tp,R2,C) Tp + log(2)*R2*C;
F.DC2 =@(Tp,T) Tp/T;
F.f =@(T) 1/T;
%====================

%////////// Values/Data /////////////////////////////////////
%............Circuit 1 (Relaxation Oscillator)
%21 = rrr, 6.7 = b,s,r, 10 = r,lb,b

C1.R1 = [47.75] * 1000; %Ohms
C1.R2 = [47.47] * 1000; %Ohms
         %100,          1000,      10000,    100000
C1.R3 = 1000*[21.75, 6.778, 2.199, .986];
C1.C  = [325e-9, 45e-9, 33.19e-9, 4.5e-9]; % Farads;


C1.DataF = [64.9, 1520, 85400, 272000 ];%Hz,kHz %Frequency used for datapoints
C1.Slew = 1.92; %Volts/us INCLUDE IN WRTIE-UP

for i = 1:length(C1.R3)
C1.TheoryF(i) = F.F1(C1.R3(i),C1.C(i));
C1.TheoryD(i) = F.y(C1.R3(i),C1.C(i));
C1.DataD(i)   = F.y(C1.R3(i),C1.C(i));
end

%plots
figure

loglog(C1.TheoryD,C1.TheoryF); hold on
loglog(C1.DataD,C1.DataF,'xr'); 
legend('Theory','Data');xlabel('RC Values');ylabel('Frequency');
hold off
title('Op-Amp Relaxation Oscillator')  
%+++++++++++++++++++++++++++++++++++++++++++++++++

%............Circuit 2 (555 Timer Chip)
C2.R1 = [.98] * 1000; %Ohms
C2.R2 = 1000 * [21.75, 6.778, 2.199, .986];%Ohms
C2.C = [325e-9, 45e-9, 33.19e-9, 4.5e-9]; %Farads

C2.R1b = 1000 * [.98, 10, 47.41, 400];
C2.R2b = [21.75] * 1000;
C2.C1b = [325.6e-6];

C2.DataF = [64.9, 1520, 85400, 272000]; %Frequency used for datapoints
C2.DataDC = [50, 63, 81, 95]; %Frequency used for datapoints

for i = 1:length(C2.R2)
C2.TheoryF(i) = F.F2(C2.R1,C2.R2(i),C2.C(i));
C2.TheoryD(i) = F.y(C2.R2(i),C2.C(i));
C2.DataD(i)   = F.y(C2.R2(i),C2.C(i));
end

for i = 1:length(C2.R1b)
   C2.TheoryDC(i) = F.DC(C2.R1b(i),C2.R2b)*100;  
end

%Plot Freq/RC
figure
subplot(2,1,1)

loglog(C2.TheoryF,C2.TheoryD);hold on
title('555 Timer Oscillator')
loglog(C2.DataF,C2.DataD,'x');
xlabel('Frequency');ylabel('RC Value')
legend('Theory','Data')
hold off

subplot(2,1,2) 
semilogx(C2.R1b,C2.TheoryDC); hold on
semilogx(C2.R1b,C2.DataDC,'x');
xlabel('R1 Values');ylabel('Duty Cycle');
legend('Theory','Data')
hold off


%+++++++++++++++++++++++++++++++++++++++++++++++++

%............Circuit 3 (Voltage controlled oscillator)
C3.R1 = [47] * 1000;
C3.R2 = [47] * 1000;
C3.C =  10e-9 ; 
C3.VValues = [.5,1.5,2.5,3.5,4.5];

C3.DataF = [2599,2000,1421,950.048655494938,518.607001671227];
C3.DataDC = [15,35,55,70,88];

for i = 1:length(C3.VValues)
    TheoryTP = F.TP(C3.R1,C3.R2,C3.C,C3.VValues(i));
    TheoryTT = F.T(TheoryTP,C3.R2,C3.C);
    C3.TheoryDC(i) = TheoryTP/TheoryTT * 100; 
    C3.TheoryF(i) = F.f(TheoryTT);
end

figure
subplot(2,1,1)
title('555 VCO')
    hold on 
    plot(C3.VValues,C3.TheoryF)
   plot(C3.VValues,C3.DataF,'x');
   legend('Theory','Data')
   xlabel('Voltage');ylabel('Oscillation Frequency');
    hold off


subplot(2,1,2)

    hold on 
    plot(C3.VValues,C3.TheoryDC)
   plot(C3.VValues,C3.DataDC,'x')
   legend('Theory','Data')
   xlabel('Voltage');ylabel('Duty Cycle');
    hold off


    
%+++++++++++++++++++++++++++++++++++++++++++++++++

%========================== ========================