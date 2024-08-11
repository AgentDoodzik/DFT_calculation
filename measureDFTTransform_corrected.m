%% amplitude * 1 - import data
close all
clear all


% File paths
filePaths = {
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__0_1Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__0_125Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__0_167Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__0_25Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__0_5Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__1Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__2Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__1\A1__4Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__0_1Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__0_125Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__0_167Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__0_25Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__0_5Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__1Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__2Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_75\A0_75__4Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__0_1Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__0_125Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__0_167Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__0_25Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__0_5Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__1Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__2Hz.csv',
    'C:\Users\macie\OneDrive\Dokumenty\Politechnika\projekt_inz\DFT_calculation\amplituda__0_5\A0_5__4Hz.csv'
};

% Read and extract columns from each file
% transposed for convenience
% first row - input (sine wave)
% second row - output
A1_sineFreq1 = table2array(readSelectedColumns(filePaths{1}))';
A1_sineFreq2 = table2array(readSelectedColumns(filePaths{2}))';
A1_sineFreq3 = table2array(readSelectedColumns(filePaths{3}))';
A1_sineFreq4 = table2array(readSelectedColumns(filePaths{4}))';
A1_sineFreq5 = table2array(readSelectedColumns(filePaths{5}))';
A1_sineFreq6 = table2array(readSelectedColumns(filePaths{6}))';
A1_sineFreq7 = table2array(readSelectedColumns(filePaths{7}))';
A1_sineFreq8 = table2array(readSelectedColumns(filePaths{8}))';

%A0_75_sineFreq1 = table2array(readSelectedColumns(filePaths{9}))'; % corrupt measurement
A0_75_sineFreq2 = table2array(readSelectedColumns(filePaths{9}))';
A0_75_sineFreq3 = table2array(readSelectedColumns(filePaths{10}))';
A0_75_sineFreq4 = table2array(readSelectedColumns(filePaths{11}))';
A0_75_sineFreq5 = table2array(readSelectedColumns(filePaths{12}))';
A0_75_sineFreq6 = table2array(readSelectedColumns(filePaths{13}))';
A0_75_sineFreq7 = table2array(readSelectedColumns(filePaths{14}))';
A0_75_sineFreq8 = table2array(readSelectedColumns(filePaths{15}))';

A0_5_sineFreq1 = table2array(readSelectedColumns(filePaths{16}))';
A0_5_sineFreq2 = table2array(readSelectedColumns(filePaths{17}))'; 
A0_5_sineFreq3 = table2array(readSelectedColumns(filePaths{18}))';
A0_5_sineFreq4 = table2array(readSelectedColumns(filePaths{19}))';
A0_5_sineFreq5 = table2array(readSelectedColumns(filePaths{20}))';
A0_5_sineFreq6 = table2array(readSelectedColumns(filePaths{21}))';
A0_5_sineFreq7 = table2array(readSelectedColumns(filePaths{22}))';
A0_5_sineFreq8 = table2array(readSelectedColumns(filePaths{23}))';

sineFrequencies = [0.0960, 0.1330, 0.1664, 0.2470, 0.5004, 1.0164, 2.0151, 4.0099]; %Hz, slightly differs from initially accepted test frequencies

matrixOfResults={A1_sineFreq1, A1_sineFreq2, A1_sineFreq3, A1_sineFreq4, ...
    A1_sineFreq5, A1_sineFreq6, A1_sineFreq7, A1_sineFreq8, ...
    A0_75_sineFreq2, A0_75_sineFreq3, A0_75_sineFreq4, A0_75_sineFreq5,...
    A0_75_sineFreq6, A0_75_sineFreq7, A0_75_sineFreq8,...
    A0_5_sineFreq1, A0_5_sineFreq2, A0_5_sineFreq3, A0_5_sineFreq4,...
    A0_5_sineFreq5, A0_5_sineFreq6, A0_5_sineFreq7, A0_5_sineFreq8 } 

AmplitudeSineFreqInfo = {'A = 1, f_{sin} = 0.0960Hz', 'A = 1, f_{sin} = 0.1330Hz', 'A = 1, f_{sin} = 0.1664Hz', 'A = 1, f_{sin} = 0.2470Hz', 'A = 1, f_{sin} = 0.5004Hz', 'A = 1, f_{sin} = 1.0164Hz', 'A = 1, f_{sin} = 2.0151Hz', 'A = 1, f_{sin} = 4.0099Hz', 'A = 0.75, f_{sin} = 0.1330Hz', 'A = 0.75, f_{sin} = 0.1664Hz', 'A = 0.75, f_{sin} = 0.2470Hz', 'A = 0.75, f_{sin} = 0.5004Hz', 'A = 0.75, f_{sin} = 1.0164Hz', 'A = 0.75, f_{sin} = 2.0151Hz', 'A = 0.75, f_{sin} = 4.0099Hz', 'A = 0.5, f_{sin} = 0.0960Hz', 'A = 0.5, f_{sin} = 0.1330Hz', 'A = 0.5, f_{sin} = 0.1664Hz', 'A = 0.5, f_{sin} = 0.2470Hz', 'A = 0.5, f_{sin} = 0.5004Hz', 'A = 0.5, f_{sin} = 1.0164Hz', 'A = 0.5, f_{sin} = 2.0151Hz', 'A = 0.5, f_{sin} = 4.0099Hz'};

% A0_75_sineFreq1 - corrupt - skipped!
%%
iter = 1;
freqSelector = 1;
columnSelector=1


MagnitudePhaseInput = {zeros(2,length(sineFrequencies)), zeros(2,length(sineFrequencies)), zeros(2,length(sineFrequencies)) }; %data from input signal 
MagnitudePhaseOutput = {zeros(2,length(sineFrequencies)), zeros(2,length(sineFrequencies)), zeros(2,length(sineFrequencies)) }; %from output signal

Gs = {zeros(1, length(sineFrequencies)), zeros(1, length(sineFrequencies)), zeros(1, length(sineFrequencies))}; % gain, 1st row - A = 1.0; 2nd - A = 0.75; 3rd A = 0.5 
fis = {zeros(1, length(sineFrequencies)), zeros(1, length(sineFrequencies)), zeros(1, length(sineFrequencies))}; % phase, rows set as above


ReWe=zeros(1,length(sineFrequencies)); %real part of input signal
ImWe=zeros(1,length(sineFrequencies)); %Imaginary part of output signal

ReWy=zeros(1,length(sineFrequencies));
ImWy=zeros(1,length(sineFrequencies));


fp= 81.4;% <- sampling frequency
fn=fp/2; %Nyquist frequency

%M = 351; %number of points used in the moving average filter 
M=351;
for iter=1:length(matrixOfResults);

matrixOfResults{iter}(2,:) =  (matrixOfResults{iter}(2,:) + ((matrixOfResults{iter}(2,:) > 50)*(-256)));
Nt = size(matrixOfResults{iter}, 2); %number of samples
%sampleNumberVector=1:Nt;
filteredSignal = zeros(1,Nt);


t=(0:Nt-1)/fp; %time vector needed for DFT

Nf = round(Nt/2+1); %sample corresponding to Nyquist frequency

%Moving average filter
for iter2=(M-1):(Nt-(M-1));
    partialSum = 0;
    for iter3=-(M-1)/2:(M-1)/2;
        partialSum = partialSum + matrixOfResults{iter}(2,iter2+iter3);
    end
    partialSum = partialSum/M;
    filteredSignal(iter2) = partialSum;
end

if(iter==9)
    freqSelector = 2; %skipping the corrupt measurement belonging to series with A = 0.75;
end

freq=sineFrequencies(freqSelector); %calculating DFT for a specific test freq.

ReWe(columnSelector) = 2/Nt * korL(matrixOfResults{iter}(1,:), cos(2*pi*freq*t)); % real part DFT 
ImWe(columnSelector) = abs(2/Nt * korL(matrixOfResults{iter}(1,:), sin(2*pi*freq*t))); %imaginary part DFT

ReWy(columnSelector) = 2/Nt * korL(filteredSignal, cos(2*pi*freq*t));
ImWy(columnSelector) = 2/Nt * korL(filteredSignal, sin(2*pi*freq*t));
    

MagnitudePhaseInput{iter}(1,columnSelector)=sqrt(1,columnSelector)^2 + ImWe(1,columnSelector)^2));
MagnitudePhaseInput{iter}(2,columnSelector)= atan2(ImWe(1, columnSelector),ReWe(1, columnSelector));

MagnitudePhaseOutput{iter}(1,columnSelector)=sqrt(ReWy(1, columnSelector)^2 + ImWy(1, columnSelector)^2);
MagnitudePhaseOutput{iter}(2,columnSelector)= atan2(ImWy(1, columnSelector),ReWy(1, columnSelector));

if (iter==8) || (iter == 16)
rowSelector1of6 = rowSelector1of6+2;
rowSelector1of3 = rowSelector1of3+1;
end

if(freqSelector ~= 8)
    freqSelector = freqSelector+1;
    columnSelector = columnSelector+1;
else
    columnSelector = 1;
    freqSelector = 1;
end

end

Gs = MagnitudePhaseOutput{:}(1,:) / MagnitudePhaseInput{:}(1,:);
fis = MagnitudePhaseOutput{:}(2,:) - MagnitudePhaseInput{:}(2,:); 


% fis = mod(fis + pi, pi);
% fis = fis - pi;

for plotIter=1:size(Gs,1)
    
    figure(plotIter);
    %title(t, AmplitudeSineFreqInfo);
    
    t = tiledlayout(2,1);
    
    nexttile
    plot(sineFrequencies,Gs(plotIter,:),'b-');
    %stem(sineFrequencies,Gs(plotIter,:));
    xlim([0 4.25]);
    title("Ch-ka Bode'ego G(f)");
    xlabel('Hz')
    ylabel('G')
    
    nexttile
    %plot(sineFrequencies,rad2deg(fis(plotIter,:)));
    stem(sineFrequencies,rad2deg(fis(plotIter,:)));
    xlim([0 4.25]);
    title("Ch-ka Bode'ego dfi(f)");
    xlabel('Hz')
    ylabel('delta fi')
    
    

end

