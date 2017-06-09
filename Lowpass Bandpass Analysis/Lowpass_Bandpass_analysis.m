%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Physics Junior Lab, 2016                           /
%     University of New Mexico                           /
%/////////////////////////////////////////////////////////

%LF411 and LM741 OP-Amps
%Lecture Notes Vout/Vin = -Z2/Z1=G(w)=N(w0/D(w)    N(w) = [(wN1
%+jw)(wN2+jw).....(wNn+jw)]
%LowPass = (-R2./R1).*(1./(1+1j.*w.*R2.*C1))
clc
clear all;
close all

%////////////////While doing Lab - Initializations//////////

%+++++++++++++++++Section 1 Lab++++++++++++++++++++++++++++
R1 = 9.91*1000       ;  %Ohms
R2 = 51.25*1000      ; %Ohms
OperatingF = 239*1000; %Function Generator Frequency
vpp = 1;
LM741_OpenLoopGain = 63000;
LF411_OpenLoopGain = 1255943;
%======================================================

%+++++++++++++++++Low Pass++++++++++++++++++++++++++++
LP_R1 = 4.506 * 1000   ; %4.7kOhms -> to Ohms 
LP_R2 = 4.552 * 1000   ; %4.7kOhms  -> to Ohms
LP_C1 = 9.49e-9   ; %10nF -> Farads

LP_DATA_Freq =[100.024000000000;146.776000000000;215.508000000000;316.277000000000;464.171000000000;681.356000000000;1000.05400000000;1467.76400000000;2154.52200000000;3162.21300000000;4641.52500000000;6812.99700000000;9999.98300000000;14678.0160000000]';
LP_DATA_AMP = [0.0520000000000000;0.0430000000000000;0.0370000000000000;0.0150000000000000;-0.0280000000000000;-0.111000000000000;-0.290000000000000;-0.637000000000000;-1.29800000000000;-2.41300000000000;-4.10200000000000;-6.40800000000000;-9.05200000000000;-12.0170000000000]';
LP_DATA_PHASE = [178.360000000000;177.606000000000;176.569000000000;175.040000000000;172.788000000000;169.439000000000;164.739000000000;158.272000000000;149.632000000000;139.659000000000;128.924000000000;119.428000000000;111.295000000000;104.997000000000]';
%======================================================

%+++++++++++++++++Bandpass++++++++++++++++++++++++++++
BP_R1 = 465;    %470 Ohms
BP_R2 = 4.552 * 1000;   %4.7kOhms -> Ohms
BP_C1 = 330e-9; %330nF -> F 
BP_C2 = 20.47e-9;   %22nF -> F

BP_DATA_Freq = (pi/2).*[100.024000000000;158.511000000000;251.271000000000;398.047000000000;630.878000000000;1000.05400000000;1584.92500000000;2511.96300000000;3981.03100000000;6309.52400000000;9999.98300000000;15848.8750000000]';  %Recorded Data
BP_DATA_AMP = [-6.60300000000000;-3.93900000000000;-2.22900000000000;-1.33400000000000;-1.00800000000000;-1.15500000000000;-1.83600000000000;-3.23000000000000;-5.40600000000000;-8.10300000000000;-10.8170000000000;-12.9750000000000]';   %Recorded Data
BP_DATA_PHASE = [-124.691000000000;-137.949000000000;-152.780000000000;-167.968000000000;177.528000000000;162.172000000000;144.309000000000;123.367000000000;100.166000000000;76.1500000000000;51.5310000000000;26.9790000000000]'; %Recorded Data
%======================================================
%//////////////////////END LAB///////////////////////////




%<<<<<<<<<<< Gain Bandwidth Product >>>>>>>>>>>>>>>>
gain = -R2/R1;
amp = 20*log10(abs(gain));
SR = -6.298; %411
SR2 = -.684 ;
closed_loop_bandwidth741 = abs(LM741_OpenLoopGain/gain);
closed_loop_bandwidth411 = abs(LF411_OpenLoopGain/gain);
string = sprintf('Part 1\n------------------------------------------');
disp(string)
string = sprintf('Gain = %g ',gain);
disp(string)
string = sprintf('Gain(dB) = %g ',amp);
disp(string)
string = sprintf('741 CLBW Peak = %g kHz',closed_loop_bandwidth741/1000);
disp(string)
string = sprintf('411 CLBW Peak = %g kHz',closed_loop_bandwidth411/1000);
disp(string)
string = sprintf('Slew Rate (LM741) = %g V/us',SR);
disp(string)
string = sprintf('Slew Rate (LF411) = %g V/us',SR2);
disp(string)
%-------------------------------------------

%<<<<<<<<<<< Low-pass Filter >>>>>>>>>>>>>>>>
%LF411

if ~isempty(LP_R1) && ~isempty(LP_R2) && ~isempty(LP_C1)
Gain_LP = -LP_R2/LP_R1;
amp = 20*log10(abs(Gain_LP));

Breakpoint_LP = 1/(LP_C1*LP_R2);
FreqEnd = 1/(LP_C1*LP_R1);
string = sprintf('\nPart 2 (Lowpass)\n------------------------------------------');
disp(string)
string = sprintf('Gain = %g',Gain_LP);
disp(string)
string = sprintf('Gain(dB) = %g',amp);
disp(string)
string = sprintf('Breakpoint = %g',Breakpoint_LP);
disp(string)
string = sprintf('Frequency Max (yaxis = 0) = %g',FreqEnd);
disp(string)
end
%-------------------------------------------


%<<<<<<<<<<< Band-pass Filter >>>>>>>>>>>>>>>>

if ~isempty(BP_R1) && ~isempty(BP_R2) && ~isempty(BP_C1) && ~isempty(BP_C2) 
passbandStart = 1/(BP_R1*BP_C1);
passbandEnd = 1/(BP_R2*BP_C2);
centerFreq = sqrt(passbandStart*passbandEnd);
passbandRange = passbandEnd-passbandStart;
Gain_BP = -BP_R2/BP_R1;
amp = 20*log10(abs(Gain_BP));
R1_Ideal = 470; R2_Ideal = 4700; %Check Values for <20% 
C1_Ideal = 330e-9; C2_Ideal = 22e-9;   %Check Values for <20% 
string = sprintf('\nPart 3 (Bandpass)\n------------------------------------------');
disp(string)
check = 100 * [BP_R1./R1_Ideal, BP_R2./R2_Ideal, BP_C1./C1_Ideal, BP_C2./C2_Ideal];
if find(check < 20)
string = sprintf('Check Values >= 20....\n %g, %g, %g, %g',check(1),check(2),check(3),check(4));
disp(string)
end
string = sprintf('Gain(dB) = %g',amp);
disp(string)
string = sprintf('Passband w<< = %g',passbandStart);
disp(string)
string = sprintf('Passband w>> = %g',passbandEnd);
disp(string)
string = sprintf('Center Frequency = %g',centerFreq);
disp(string)
end
%-------------------------------------------    

finish = 1;

clear('ans','D','band','Breakpoint_LP','C1_Ideal','C2_Ideal','R1_Ideal','R2_Ideal'...
,'OperatingF','R1','R2','SR','string','vpp','X_c','w','X_c1','X_c2','Z_f','Z_in',...
'LF411_OpenLoopGain','LM741_OpenLoopGain','Gain_LP',...
'Gain_BP','gain','FreqEnd','fr','check','closed_loop_bandwidth411','closed_loop_bandwidth741',...
'centerFreq','C1','C2','amp','OperatingF_LP','SR2')
%% <<<<< Write-Up >>>>>>>>
if finish >= 1
%****Initialization (All Anonymous Functions Included)****
f = logspace(0,5,10000); %
w = 2*pi*f;
%Lowpass
Z_f_LP =@(R1,R2) R2./R1;               
Z_in_LP =@(w,R2,C1) (1+1j.*w.*R2.*C1); 
%Bandpass
Z_f_BP =@(w,R2,C1) (1j.*w.*R2.*C1);    
Z_in_BP =@(w,R1,C1,R2,C2) ((1j.*w.*R1.*C1 + 1).*(1j.*w.*C2.*R2 + 1)); 
%Z_f_BP =@(w,R,C) (1./(1./R + 1j.*w.*C));
%Z_in_BP =@(w,R,C) (R + 1./(1j.*w.*C));


amp =@(G) 20.*log10(abs(G));
phaseA =@(G) (180/pi).*atan(imag(G)./real(G));

PassbandStart =@(R1,C1) 1/(2*pi*R1*C1); %First passband *Gives Frequency
PassbandEnd =@(R2,C2) 1/(2*pi*R2*C2);   %Second passband *Gives Frequency
PassbandCenter =@(P1,P2) sqrt(P1*P2);       %Center of Passband *Gives Frequency
%***********************

%//////////////////////////  Low-pass /////////////////////////////////////

LP_G = -Z_f_LP(LP_R1,LP_R2) ./ Z_in_LP(w,LP_R2,LP_C1) ;
LP_Amp = amp(LP_G);
LP_Phase = phaseA(LP_G);
LP_Band = PassbandStart(LP_R2,LP_C1);
LP_Amp_Max = max(LP_Amp);
LP_Amp_Min = min(LP_Amp);

figure('units','normalized','position',[.1 .1 .4 .8])
%PLOTS LP AMPLITUDE -----------
subplot(2,1,1)
semilogx(f,LP_Amp);hold on
title('Lowpass Filter - Amplitude')
scatter(LP_DATA_Freq, LP_DATA_AMP,'rx')
line([LP_Band,LP_Band],[min(LP_Amp),-3],'Color','m')
line([1,LP_Band],[-3,-3],'Color','m')
text(LP_Band,-3,'      \leftarrow Band @ -3dB, 3684.3Hz')
hold off
xlabel('Frequency(Hz)');ylabel('Amplitude(dB)');legend('Theory','Data')


%PLOTS LP PHASE -----------
subplot(2,1,2)
semilogx(f,LP_Phase,'b');hold on
title('Lowpass Filter - Phase')
scatter(LP_DATA_Freq, LP_DATA_PHASE-180,'r*');hold off
xlabel('Frequency(Hz)');ylabel('Phase(\phi)');legend('Theory','Data')

%-----------------------------

%//////////////////////////  Band-pass ////////////////////////////////////

BP_G = -Z_f_BP(w,BP_R2,BP_C1) ./ Z_in_BP(w,BP_R1,BP_C1,BP_R2,BP_C2) ; 
BP_Amp = amp(BP_G);
BP_Phase = phaseA(BP_G);
BP_Amp_Max = max(BP_Amp);
BP_Amp_Min = min(BP_Amp);
BP_Pass1 = PassbandStart(BP_R1,BP_C1);
BP_Pass2 = PassbandEnd(BP_R2,BP_C2);
BP_Center = PassbandCenter(BP_Pass1,BP_Pass2);



figure('units','normalized','position',[.1 .1 .4 .8])
%PLOTS BP AMPLITUDE -----------

subplot(2,1,1) 
semilogx(f,BP_Amp,'b');
hold on
title('Bandpass Filter - Amplitude')
scatter(BP_DATA_Freq, BP_DATA_AMP,'rx')
plot([min(f),BP_Pass1],[BP_Amp(2),BP_Amp_Max],'Color','m') %Slope from left to right
line([BP_Pass2,max(f)],[BP_Amp_Max,BP_Amp(end)],'Color','m') %Slope from right to left
line([BP_Pass1,BP_Pass1],[BP_Amp_Min,BP_Amp_Max],'Color','m') %First Band
line([BP_Pass2,BP_Pass2],[BP_Amp_Min,BP_Amp_Max],'Color','m') %Second Band
line([BP_Pass1,BP_Pass2],[BP_Amp_Max,BP_Amp_Max],'Color','m') %Horizontal band (Range from band1 to band2)
line([BP_Center,BP_Center],[BP_Amp_Max,BP_Amp_Max],'Color','m','LineStyle',':') %Center of bands Line
xlabel('Frequency(Hz)');ylabel('Amplitude(dB)');legend('Theory','Data','Approximate')
hold off
%-----------------------------

%PLOTS BP PHASE -----------
subplot(2,1,2)
semilogx(f,BP_Phase,'b');
hold on
title('Bandpass Filter - Phase')
scatter(BP_DATA_Freq, unwrap(BP_DATA_PHASE)-90,'r*');
hold off
xlabel('Frequency(Hz)');ylabel('Phase(\phi)');legend('Theory','Data')
hold off
%-----------------------------


%//////////////////////////  Part 2  ////////////////////////////////////
Fc0 = 1000;   % Hz (Given)
Fc1 = 10000;  % Hz (Given)
Fc2 = 100000; % Hz (Given)
Tau1 = 1000;   % Hz (Given)
Tau2 = 10000;  % Hz (Given)
Tau3 = 100000; % Hz (Given)
f = logspace(0,10,100000);
w = 2*pi*f;
PB_R1 = 1000; %10kOhms Chosen
PB_R2 = 1000; %10kOhms Chosen

%Calculated Values C1,C2, by choosing R1,R2 
PB_C1 = 1/(2*pi*Fc1*PB_R1);
PB_C2 = 1/(2*pi*Fc2*PB_R2);
%-----------------------------------------

%Calculates Gain Peak and Bandwidth


PB_G = -Z_f_BP(w,PB_R2,PB_C1)./ Z_in_BP(w,PB_R1,PB_C1,PB_R2,PB_C2);
PB_Amp = amp(PB_G);
PB_Phase = 180+phaseA(PB_G);
PB_Gain = PB_R2/PB_R1;
PB_Gain2 = 20*log10(PB_R2/PB_R1);


PB_Amp_Max = max(PB_Amp);
PB_Amp_Min = min(PB_Amp);
PB_Band1 = PassbandStart(PB_R1,PB_C1);
PB_Band2 = PassbandEnd(PB_R1,PB_C1);
PB_Center = PassbandCenter(100000,10000);
PB_FirstLine = 1/(2*pi*PB_R2*PB_C1);
figure('units','normalized','position',[.1 .1 .4 .8])

%PART 2 AMPLITUDE -----------
subplot(2,1,1)
semilogx(f,PB_Amp,'b');hold on
title('Experimental Bandpass Filter - Amplitude')
plot([min(f),Fc1],[PB_Amp(1),PB_Amp_Max],'Color','r') %Slope from left to right
line([Fc2,max(f)],[PB_Amp_Max,PB_Amp(end)],'Color','r') %Slope from right to left
line([Fc1,Fc1],[PB_Amp_Min,PB_Amp_Max],'Color','r') %First Band
line([Fc2,Fc2],[PB_Amp_Min,PB_Amp_Max],'Color','r') %Second Band
line([Fc1,Fc2],[PB_Amp_Max,PB_Amp_Max],'Color','r') %Horizontal band (Range from band1 to band2)
line([PB_Center,PB_Center],[PB_Amp_Min,PB_Amp_Max],'Color','r','LineStyle',':') %Center of bands Line
line([Fc0,Fc0],[-20,PB_Amp_Min],'Color','r') 
text(Fc0,-30,'\leftarrow -20dB @ 1000Hz')
text(Fc1,-40,'\leftarrow Band @ 10000Hz')
text(Fc2,-50,'\leftarrow Band @ 100000Hz')
text(PB_Center,-60,'\leftarrow Amplitude Max @ -1dB')

xlabel('Frequency(Hz)');ylabel('Amplitude(dB)');legend('Theory','Approximate')
hold off
%-----------------------------

%PART 2 PHASE -----------
subplot(2,1,2)
semilogx(f,PB_Phase,'b'); hold on
title('Experimental Bandpass Filter - Phase')
xlabel('Frequency(Hz)');ylabel('Phase(\phi)');legend('Theory')
hold off
end
%-----------------------------
%% Latex Stuff For Analysis
% Section 1 - Finding the gain-bandwidth
% Gain $= -5.17154$
% Gain(dB) $= 14.2724$ 
% LM471 Gain-Bandwidth Product $\leftarrow 63000$Hz or $.063$MHz
% LM471 Slew Rate $\leftarrow |-.684|=.684$ Volts/$\mu$s
% LF411 Gain-Bandwidth Product $\leftarrow 1255943$Hz or $1.26$MHz
% LF411 Slew Rate $\leftarrow |-6.2981|=6.3$ Volts/$\mu$s

%Lowpass Filter
%Resistor #1 $= 4506 \Ohms$
%Resistor #2 $= 4552 \Ohms$
%Capacitor   $= 9.49 \times 10^{-9)$F
%Gain $= -1.01021$
%Gain(dB) $= 0.0882213$

%Bandpass Filter
%Resistor #1 $= 465 \Ohms$
%Resistor #2 $= 4552 \Ohms$
%Capacitor #1   $= 3.3 \times 10^{-7)$F
%Capacitor #2   $= 2.047 \times 10^{-8)$F
%Gain $= 9.7892$
%Gain(dB) $= 15.6934$
%Passband #1 $= \frac{1}{465\Ohms 3.3e-07F}$ 
%Passband #2 $= \frac{1}{4552\Ohms 2.047e-08F}$

%Experimental Bandpass Filter
%Chose $R_1 = R_2 = 10$k$\Omega$
%$C_1 = 1.592e-08$F
%$C_2 = 1.592e-09$F

%GENERAL Functions Used
%Lowpass
% \[Z_f = \frac{R_2}{R_1}\\ Z_in = 1+j\omega R_2C_1)\\ G(\omega) = -\frac{R2}{R1}\frac{1}{\j\omega R_2C_1 +1}\]
%Bandpass
%\[Z_f = R_2||\frac{1}{j\omega C_2}\\ Z_in = R_1 + \frac{1}{j\omega C_1}\\ G(\omega) = -\frac{j\omega R_2C_1}{(j\omega R_1C_1+1)(j\omegaC_2R_2 +1)} 


