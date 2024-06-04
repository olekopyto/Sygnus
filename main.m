%{
A pulse signal processing utility,
 given a CSV file with two columns (time [s], amplitude [v]),
 gives: Signal frequency and phase spectrum,  pulse rise, fall times,
 pulse lenght, pulse energy, area under pulse.
 It prints out the first 40 harmonics of the signal to a .csv file and
 also constructs a .asc file for implementing this model in LTSpice.
 This file contains 40 voltage sources connected together in series,
 aproximating the input signal. Originaly made for a CERN project.
 Aleksander Kopyto, AGH University of Kraków, April 2024.
%}

function [frequencies, P1_dB, phase, P2]= fftLog(signal, time)
    % N to długość sygnału odczytanego
    N = length(signal);
    
    % obliczenie transformaty fouriera
    Y = fft(signal,N);
    P2 = abs(Y/N); 
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);  
    Fs = 1 / mean(diff(time));  
    f = (Fs * (0:(N/2)) / N) / 1e6;  % Oś częstotliwości w MHz
    % Spektrum częstotliwosći wyskalowane w dB
    
    phaseDoubleSide = angle(Y);
    phaseRad = phaseDoubleSide(1:floor(N/2)+1);
    phase = unwrap(phaseRad)*(180/pi); %unwrap(phaseRad)*
    P1_dB = 20 * log10(P1);
    P1_dB(isinf(P1_dB)) = -160;  % Uzyskane w transformacie wartośći -Inf zrzuca do -160 dB
    
    frequencies = f;
    
end

function bandwidth = SpectrumAccuracy(spectrum, frequencies, desiredAccuracy)
    spectrumSum = 0;

     for i = 1:length(frequencies)
        spectrumSum = spectrumSum +  spectrum(i);   
    end

    spectrumSubSum = 0;
    bandwidth = 0;
    display(spectrumSum);
    for i = 1:length(frequencies)
        spectrumSubSum = spectrumSubSum + spectrum(i);
        
        fprintf('Frequency: %f, SpectrumSubSum: %f\n', frequencies(i), spectrumSubSum);
        if spectrumSubSum >= spectrumSum * (1-desiredAccuracy)
            bandwidth = frequencies(i);
            break;
        end
    end
end





%program pyta użytkownika o nazwę pliku

harmonicsN = 400; %liczba 
prompt = {'Proszę podać nazwę pliku csv, w folderze skryptu, do odczytu.'};
dlgtitle = 'Nazwa wymagana.';
dims = [1 35];
defaultin = "data.csv";
%defaultin = "HV1200.csv";
filename=inputdlg(prompt, dlgtitle, dims, defaultin);
% Odczyt danych z pliku csv, DOMYŚLNIE data.csv w lokalizacji skryptu.
dataTable = readtable(filename{1});
%czyści tablicę z błędnych rekordów
%w jedym pliku pojawił się rekord z czasem =0, reszta czasów z reguły
% > 0,15 us, usuwamy je, żeby nie mieliśmy błędnie rysowanych wykresów i
% bezpodstawnych przerw w impulsie, które prowadzą do błędnej normalizacji. 
dataTable = dataTable(dataTable.Var1 > 0, :);

head(dataTable);

% pobranie danych z tablicy do wektorów matlabowych
time = dataTable.Var1;
signalOrg = dataTable.Var2;

%inicjalizacja
maxAmp=0;
maxIte = 0;
maxTime = 0;
currAmp = 0;
edgeAmp = 0;
edgeRise1 = 0;
edgeRise2 = 0;
edgeFall1 = 0;
edgeFall2 = 0;
offsetDC = 0;
%%wyliczanie przedimpulsowego offsetu DC
for i = 1:floor(length(time)/10)
    offsetDC = offsetDC + signalOrg(i);
end
offsetDC = offsetDC/floor(length(time)/10);
signalNoOffset = signalOrg - offsetDC;
%%poszukiwanie amplitudy maksymalnej
for i = 1:length(time)
    currAmp = signalNoOffset(i);
    if(abs(currAmp)>abs(maxAmp))
        maxAmp = currAmp;
        maxIte = i;
        maxTime = time(i);
    end
end
disp(maxAmp);
%obliczanie czasu trwania i energii impulsu
%czas trwania definiujemy jako czas kiedy amplituda impulsu
%jest wyższa niż pół amplitud maksymalnej
halfMaxAmp = abs(0.5*maxAmp);
disp(halfMaxAmp);
pulseWidthI = find(abs(signalNoOffset) > halfMaxAmp);
pulseWidthTime = time(pulseWidthI(end)) - time(pulseWidthI(1));
%zakładamy, że dostarczony sygnał jest napięciem.
%tutaj energia jest całką kwadratu amplitudy po czasie, podzieloną przez opór.
%Przedstawiam dla oporów 1 Ohm i 50 Ohm.
energy1Ohm = trapz(time(pulseWidthI(1):pulseWidthI(end)), signalNoOffset(pulseWidthI(1):pulseWidthI(end)).^2);
energy50Ohm = energy1Ohm/50;
%pole pod pulsem
area = trapz(time(pulseWidthI(1):pulseWidthI(end)), signalNoOffset(pulseWidthI(1):pulseWidthI(end)));

%normalizacja amplitudy
signal = signalNoOffset/maxAmp;

%obliczanie czasu wzrastania/opadania
for i = 4:length(time)
    if(signal(i)>0.1 && edgeRise1==0)
        edgeRise1 = i;
    end
    if(signal(i)>0.9 & edgeRise2==0)
        edgeRise2 = i;
        break;
    end
end

for i = length(time):-1:4
    if(signal(i)>0.1 && edgeFall2==0)
        edgeFall2 = i;
    end
    if(signal(i)>0.9 & edgeFall1==0)
        edgeFall1 = i;
        break;
    end
end

timeRiseTime = time(edgeRise2) - time(edgeRise1);
timeFallTime = time(edgeFall2) - time(edgeFall1);

disp("Edges");
disp(edgeRise1);
disp(edgeRise2);
disp(edgeFall1);
disp(edgeFall2);

pulseLength = time(edgeRise2) - time(edgeRise1);
timePulse = time(edgeRise1:edgeFall2);
signalPulse = signal(edgeRise1:edgeFall2);
disp(pulseLength);

%%obliczanie wycinka pliku, w którym mamy sygnał

windowStart = 1;
windowEnd = 2*maxIte;
disp(2*maxIte)
disp(length(signal))
if(2*maxIte>length(signal))
    windowEnd = length(signal);
    windowStart = 2*maxIte - length(signal);
    disp("impule nearing data end");
end
disp("impule window set.")

disp("Window:")
disp(windowStart)
disp(windowEnd)
towindowPulse = signal(windowStart:windowEnd);
towindowTime = time(windowStart:windowEnd);

disp("window lenth and timespan:")
disp(length(towindowPulse))
disp((towindowTime(4)))
disp(towindowTime(length(towindowPulse)));
disp(towindowTime(length(towindowPulse))-(towindowTime(4)));

%%okno blackmana
blackmanPulse = towindowPulse .* blackman(length(towindowPulse));

%%normalizacja amplitudy po oknie blackmana (chcemy, aby 
blackmanNormalisationFactor = length(towindowPulse)/sum(blackman(length(towindowPulse)));

blackmanPulse = blackmanPulse*blackmanNormalisationFactor;
%blackmanPulse = blackmanPulse+offsetDC;
%%transformata FFT
[f, P1_dB, phase] = fftLog(signal, time);
[blackmanf, blackmanSpectrum, blackmanPhase, blackmanSpectrumLinear] = fftLog(blackmanPulse, towindowTime);

%Wyznaczenie pasma -3dB

blackmanfMax = max(blackmanSpectrumLinear);
blackmanf3dB = blackmanfMax/sqrt(2);
blackmanf3dBLog = max(blackmanSpectrum) - 3;
indices = find(blackmanSpectrum>=blackmanf3dBLog);
blackman3dBMin=blackmanf(indices(1));
blackman3dBMax=blackmanf(indices(end));



%wyznaczenie pasm dokładności
bit10 = 0.001;
bit12 = 0.00025;
bit14 = 0.00006;

bandwidth10 = SpectrumAccuracy(blackmanSpectrumLinear, blackmanf, bit10);
bandwidth12 = SpectrumAccuracy(blackmanSpectrumLinear, blackmanf, bit12);
bandwidth14 = SpectrumAccuracy(blackmanSpectrumLinear, blackmanf, bit14);
display(bandwidth10);
display(bandwidth12);
display(bandwidth14);

% multiplot.........................................................
figure;
sgtitle("Analiza "+filename{1});

%{ 
wykresy funkcji bez okna blackmana

% Spektrum w dB
subplot(3, 2, 1);
semilogx(f, P1_dB, 'o-');
title('Spectrum sygnału [dB]');
xlabel('Częstotliwość [MHz]');
ylabel('Natężenie [dB]');

% Wykres Faz
sgtitle('Fazy z transformaty orginalnego pliku impulsu');

subplot(3, 2, 2);
plot(f, phase);
title('Sygnał');
xlabel('Częstotliwość [MHz]');
ylabel('Faza [deg]');

%}

% Spektrum w dB
subplot(2, 2, 1);

semilogx(blackmanf, blackmanSpectrum);
grid minor;
title('Spectrum sygnału [dB]');
xlabel('Częstotliwość [MHz]');
ylabel('Natężenie [dB]');
title('Spektrum częst. z transformaty FFT');
x3dbY = yline(blackmanf3dBLog, 'r--', blackmanf3dBLog+ " dBV, -3dB cutoff");
x3dbY.LabelVerticalAlignment = 'bottom';
x3dbY.LabelHorizontalAlignment = 'center';
x3dbF = xline(blackman3dBMax,'r--',blackman3dBMax+ " MHz");
x3dbF.LabelVerticalAlignment = 'middle';
x10 = xline(bandwidth10,'b--',bandwidth10 + " MHz, 10bit");
x10.LabelHorizontalAlignment = 'left';
x10.LabelVerticalAlignment = 'middle';
x12 = xline(bandwidth12,'b--',bandwidth12+ " MHz, 12bit");
x12.LabelHorizontalAlignment = 'left';
x12.LabelVerticalAlignment = 'middle';
x14 = xline(bandwidth14,'b--',bandwidth14+ " MHz, 14bit");
x14.LabelHorizontalAlignment = 'right';
x14.LabelVerticalAlignment = 'middle';
% Wykres sygnału
subplot(2, 2, 2);
plot(towindowTime, blackmanPulse);
grid minor;
title('Sygnał po normalizacji');
xlabel('Czas [s]');
ylabel('Natężenie [V]');

% Wykres Faz

subplot(2, 2, 3);
semilogx(blackmanf, blackmanPhase);
grid minor;
title('Spektrum faz z transformaty FFT');
xlabel('Częstotliwość [MHz]');
ylabel('Faza [deg]');

% Wykres orginalnego sygnału

subplot(2, 2, 4);

plot(time, signalNoOffset);
grid minor;
title('Sygnał Orginalny');
xlabel('Czas [s]');
ylabel('Natężenie');
xImpBeg = xline(time(pulseWidthI(1)),'r--',"Pocz. całkowania");
xImpBeg.LabelHorizontalAlignment="left";
xImpBeg.LabelVerticalAlignment="bottom";
xImpEnd = xline(time(pulseWidthI(end)),'r--',"Koniec całkowania");
xImpEnd.LabelVerticalAlignment="bottom";
xImpEnd.LabelHorizontalAlignment="right";

%%zapis wyników analizy do pliku....................................



outFile = fopen(filename{1}+"_out.txt","w+");
fprintf(outFile,"%s\n",datetime("now"));
fprintf(outFile,"%s %g\n","Czas wzrastania [s]:", timeRiseTime);
fprintf(outFile,"%s %g\n","Czas opadania [s]:", timeFallTime);
fprintf(outFile,"%s %g\n","Czas trwania impulsu [s]:", pulseWidthTime);
fprintf(outFile,"%s %g\n","Maksymalna amplituda impulsu [V]:", maxAmp);
fprintf(outFile,"%s %g\n","Chwila maksimum sygnału [s]:", maxTime);
fprintf(outFile,"%s %g %g\n","Energia impulsu dla 1 Ohm i 50 Ohm impedancji odbiornika impulsu [J]:", energy1Ohm, energy50Ohm);
fprintf(outFile,"%s %g\n","Pole powierzchni impulsu [V*s]:", area);
fprintf(outFile,"%s %g \n","Pasmo -3dB [MHz]: ", blackman3dBMax);
fprintf(outFile,"%s %g, %g, %g \n","Pasma dokładności (10, 12, 14 bit) [MHz]: ", bandwidth10, bandwidth12, bandwidth14);
fprintf(outFile,"Tableka harmonicznych: \n");
fprintf(outFile,"%s\n","Częstotliwość [Hz], amplituda [dBV], amplituda [V], faza [°]");
for i=1:harmonicsN
    fprintf(outFile,"%gM,%g,%g\n",blackmanf(i),blackmanSpectrum(i),mod(blackmanPhase(i),360));
end
fprintf(outFile,"%s\n"," Bez okna blackmana:");
fprintf(outFile,"%s\n"," freq, amplitude [dB], phase []");
for i=1:harmonicsN
    fprintf(outFile,"%gM,%g,%g\n",f(i),P1_dB(i),mod(phase(i),360));
end

fclose(outFile);

%dodatkowo, zapis tabeli harmonicznych do pliku csv, dla łatwiejszego
%przetwarzania, np. w Excelu czy LibreOffice Calc.
csvFile = fopen(filename{1}+"_harmonicTab.csv","w");
if(csvFile == -1)
    error("Proszę zamknąć program korzystający z pliku "+filename{1}+"_harmonicTab.csv");
end
fprintf(csvFile,"%s\n","freq [Hz], amplitude [dBV], linear amplitude [V], phase [°]");
for i=1:harmonicsN
    fprintf(csvFile,"%f,%f,%f,%f\n",10E5*blackmanf(i),blackmanSpectrum(i),blackmanSpectrumLinear(i),mod(blackmanPhase(i),360));
end
fclose(csvFile);

libFile = fopen("AliceSignal.sub","w"); %filename{1}+".sub"
if(libFile == -1)
    error("Proszę zamknąć program korzystający z pliku "+filename{1}+".asc");
end

header = ".subckt AliceSignal p n";
fprintf(libFile,"%s\n",header);

for i=2:harmonicsN
    if(blackmanSpectrumLinear(i)<0.000001)
        continue;
    end
    fprintf(libFile,"V%d p_aux%d 0 SIN(0 %f %fMeg 0 0 %f)\n",i,i,blackmanSpectrumLinear(i),blackmanf(i),blackmanPhase(i)+90);

end

fprintf(libFile,"B1 p n V = %f",blackmanSpectrumLinear(i));

for i=2:harmonicsN
    if(blackmanSpectrumLinear(i)<0.000001)
        continue;
    end
    fprintf(libFile," + V(p_aux%d,0)",i);

end
fclose(libFile);
%generacja pliku symulacji do LTSpice................................

%{

ascFile = fopen(filename{1}+"a.asc","w+");

if(ascFile == -1)
    error("Proszę zamknąć program korzystający z pliku "+filename{1}+".asc");
end

header = [
    "Version 4",
    "SHEET 1 944 680",
    ];

wire1 = "WIRE 128 160 ";
wire2 = " 160";

flag1 = "FLAG 128 ";
flag2 = " 0";

symbol1 = "SYMBOL voltage 128";
symbol2 = [ " R0",
"WINDOW 123 0 0 Left 0",
"WINDOW 39 24 124 Left 2",
"SYMATTR InstName V"
    ];
symbol2 = join(symbol2,newline);
symbol3 = newline+"SYMATTR Value SINE(0";
symbol3b = newline+"SYMATTR Value ";
symbol4 = newline+"SYMATTR SpiceLine Rser=0"; %%ok

%zapis do pliku symulacji do LTSpice................................
fprintf(ascFile,"%s\n",header);

fprintf(ascFile,"%s%d%s\n",flag1,(harmonicsN+1)*80+16,flag2);
fprintf(ascFile,"%s %d %s%d %s %f %s\n\n",symbol1,80,symbol2,1,symbol3b,blackmanSpectrumLinear(i),symbol4);
for i=2:harmonicsN
    if(blackmanSpectrumLinear(i)<0.000001)
        continue;
    end
    fprintf(ascFile,"%s %d %s%d %s %f %fMeg 0 0 %f) %s\n\n",symbol1,i*80,symbol2,i,symbol3,blackmanSpectrumLinear(i),blackmanf(i),blackmanPhase(i)+90,symbol4);
end

fclose(ascFile);

%}



disp( ...
    "Analiza skończona, wyniki zapisano do " ...
    +filename{1} ...
    +"_out.txt" ...
    +newline ...
    +"tablicę harmonicznych do: " ...
    +filename{1} ...
    +"_harmonicTab.csv" ...
    +newline ...
    +"a model impulsu dla LTSpice do: " ...
    +filename{1} ...
    +".asc" ...
);

