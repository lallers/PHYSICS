%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Physics Junior Lab, 2016                           /
%     University of New Mexico                           /
%/////////////////////////////////////////////////////////

clear all; clc; close all;
%<<Initializations Frequency/Calculated Omega>>
freq = 1:1:100000;
omega = 2*pi*freq;

%<<Circuit Components>>
highR = 1000; % Highpass Resistor Ohms
highC = 0.00000001; % Highpass Capacitor Farads
R_rcl = 471; % Ohms
C_rcl = 2.23e-10; % Farads
L_rcl = .096; % Hertz
%$$$$$$$$$$$$$$$$$$$$$$$$

%!!!!!!!!!!!!Need to add to graph still?!!!!!!!!!!!!!!!!!!!!!
w0_highRC = 1/(highR*highC);%Omega of peak frequency 1/sqrt(LC)
f0_highRC = w0_highRC/(2*pi); %Frequency Peak
w0_rcl = 1/(L_rcl*C_rcl)^(1/2); %Omega of peak frequency 1/sqrt(LC)
f0_rcl = w0_rcl/(2*pi); %Frequency Peak
%!!!!!!!!!!!!Need to add to graph still?!!!!!!!!!!!!!!!!!!!!!


%---------------------------------data Highpass RC----------------
highRCData_freq = [0.931,1.676,2.421,3.912,6.333,10.058,15.832,25.146,39.861...
    63.144,100.024,158.511,251.271,398.047,630.878,1000.054];
highRCData_amp = (-1)*[11.970,7.592,5.235,2.926,1.460,.702,.354,.179,.088...
    .040,.022,.014,.002,.006,.002,.007];
highRCData_phase = (-1)*[73,63,54,41,29,20,13,8,5,3,2,1,1,0,0,0];
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%---------------------------------data RCL Series--------
SeriesData_f = (1/10)*[9999.983,12589.246,15848.875,19952.655,25118.887,31622.685...
39810.687,50118.752,63095.801,79432.875,100000.016,125892.460,158489.309...
199526.176,251188.688,316227.786,398107.246,501187.146,630957.261,794328.190,999999.978];
SeriesData_amp = (-1)*[38.305,36.375,33.956,31.376,28.391,24.660,18.810,6.625,17.543,23.715...
    27.631,30.565,33.182,35.501,37.638,39.704,42.447,44.422,46.243,48.211,50.493];
SeriesData_phase = [88,86,86,86,84,82,75,9,-73,-81,-84,-84,-85,-86,-87,-88,-88,-87.6,-86.993,-86.483,-86.863];
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%--------------Anonymous Functions-------------------------
%<<Highpass RC Function>>
highRC_G =@(omega,R,C) R ./ (R + (1 ./ (1j.*omega.*C)));

%<<Parallel Function>>
ParallelRCL_G =@(omega,R,C,L)(1j.*omega.*L) ./ (R.*(1-(omega.^2.*L.*C)) + 1j.*omega.*L);

%<<Series Functions>>
SeriesRCL_G =@(omega,R,C,L) R./(R + (1j.*omega.*L) + (1./(1j.*omega.*C)));

%<<Convert Amplitude to Log/Phase Function>>
AmplitudeFunc =@(G) 20.*log10(abs(G));
PhaseFunc =@(G) radtodeg(atan(imag(G)./real(G)));
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%$$$$$$$$$$$$$$$$$$$$$$$$ HIGHPASS CIRCUIT $$$$$$$$$$$$$$$$$$$$
%---Theory----
highTheoryG = highRC_G(omega,highR,highC);
highAmp = AmplitudeFunc(highTheoryG);
highPhase = radtodeg(PhaseFunc(highTheoryG));

%++++++++++++++
%Plots
figure('units','normalized','position',[.1 .1 .4 .8])
subplot(2,1,1);
hold on
plot(freq,highAmp);title('Highpass Circuit (RC)');xlabel('Frequency (Hz)');ylabel('Gain (dB)');
%<<My Data>>
plot(highRCData_freq,highRCData_amp,'r*')
%-----------
%plot([highQ highQ],[-3,-3]) Experiment
hold off
subplot(2,1,2);
hold on
plot(freq,highPhase);xlabel('Frequency(Hz)');ylabel('Phase (\phi)');
%<<My Data>>
plot(highRCData_freq,highRCData_phase,'rx')
%-----------
hold off
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SERIES RCL $$$$$$$$$$$$$$$$$$$$$$$$
%Theory
Series_G = SeriesRCL_G(omega,R_rcl,C_rcl,L_rcl);
Series_Amp = AmplitudeFunc(Series_G);
Series_Phase = radtodeg(PhaseFunc(Series_G));
%------
%Plots
figure('units','normalized','position',[.1 .1 .4 .8])
subplot(2,1,1);
hold on
plot(freq,Series_Amp);title('Series Circuit (RCL)');xlabel('Frequency (Hz)');ylabel('Gain (dB)');
%<<My Data>>
plot(SeriesData_f,SeriesData_amp,'r*')
%-----------
hold off
subplot(2,1,2);
hold on
plot(freq,Series_Phase);xlabel('Frequency(Hz)');ylabel('Phase (\phi)');
%<<My Data>>
plot(SeriesData_f,SeriesData_phase,'rx')
%-----------
hold off
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$ PARALLEL RCL $$$$$$$$$$$$$$$$$$$$$$$$
%Theory
Parallel_Q = 1/((L_rcl*C_rcl)^(1/2)); %Not Graphed, No Data
Parallel_Amp = AmplitudeFunc(Parallel_Q); %Not Graphed, No Data
Parallel_Phase = radtodeg(PhaseFunc(Parallel_Q)); %Not Graphed, No Data
%theory_Q = 

figure('units','normalized','position',[.1 .1 .4 .8])
R = R_rcl; C = C_rcl; L = L_rcl; G1 = tf([1/(R*C) 0],[1 1/(R*C) 1/(L*C)]);
R2 = 1000+R_rcl; G2 = tf([1/(R2*C) 0],[1 1/(R2*C) 1/(L*C)]);
R3 = 1000-R_rcl; G3 = tf([1/(R3*C) 0],[1 1/(R3*C) 1/(L*C)]);
R4 = 10000; G4 = tf([1/(R4*C) 0],[1 1/(R4*C) 1/(L*C)]);
R5 = 100000; G5 = tf([1/(R5*C) 0],[1 1/(R5*C) 1/(L*C)]);
%hold on
bode(G1,'b',G2,'g',G3,'r',G4,'y',G5,'m');title('Parallel Circuit (RCL)')
%subplot(1,1,1); plot(cQx,cQy,'rx')
%hold off
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

