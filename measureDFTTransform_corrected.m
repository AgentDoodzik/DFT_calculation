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
A1__0_1Hz = table2array(readSelectedColumns(filePaths{1}))';
A1__0_125Hz = table2array(readSelectedColumns(filePaths{2}))';
A1__0_167Hz = table2array(readSelectedColumns(filePaths{3}))';
A1__0_25Hz = table2array(readSelectedColumns(filePaths{4}))';
A1__0_5Hz = table2array(readSelectedColumns(filePaths{5}))';
A1__1Hz = table2array(readSelectedColumns(filePaths{6}))';
A1__2Hz = table2array(readSelectedColumns(filePaths{7}))';
A1__4Hz = table2array(readSelectedColumns(filePaths{8}))';

% A0_75__0_1Hz = table2array(readSelectedColumns(filePaths{9}))';
% A0_75__0_125Hz = table2array(readSelectedColumns(filePaths{10}))';
% A0_75__0_167Hz = table2array(readSelectedColumns(filePaths{11}))';
% A0_75__0_25Hz = table2array(readSelectedColumns(filePaths{12}))';
% A0_75__0_5Hz = table2array(readSelectedColumns(filePaths{13}))';
% A0_75__1Hz = table2array(readSelectedColumns(filePaths{14}))';
% A0_75__2Hz = table2array(readSelectedColumns(filePaths{15}))';
% A0_75__4Hz = table2array(readSelectedColumns(filePaths{16}))';
% 
% A0_5__0_1Hz = table2array(readSelectedColumns(filePaths{17}))';
% A0_5__0_125Hz = table2array(readSelectedColumns(filePaths{18}))';
% A0_5__0_167Hz = table2array(readSelectedColumns(filePaths{19}))';
% A0_5__0_25Hz = table2array(readSelectedColumns(filePaths{20}))';
% A0_5__0_5Hz = table2array(readSelectedColumns(filePaths{21}))';
% A0_5__1Hz = table2array(readSelectedColumns(filePaths{22}))';
% A0_5__2Hz = table2array(readSelectedColumns(filePaths{23}))';
% A0_5__4Hz = table2array(readSelectedColumns(filePaths{24}))';

sineFrequencies = [0.0960, 0.1330, 0.1664, 0.2470, 0.5004, 1.0164, 2.0151, 4.0099]; %Hz, slightly differs from initially accepted test frequencies

matrixOfResults={A1__0_1Hz, A1__0_125Hz, A1__0_167Hz, A1__0_25Hz, A1__0_5Hz, A1__1Hz, A1__2Hz, A1__4Hz} %, A0_75__0_1Hz, A0_75__0_125Hz, A0_75__0_167Hz, A0_75__0_25Hz, A0_75__0_5Hz, A0_75__1Hz, A0_75__2Hz, A0_75__4Hz, A0_5__0_1Hz, A0_5__0_125Hz, A0_5__0_167Hz, A0_5__0_25Hz, A0_5__0_5Hz, A0_5__1Hz, A0_5__2Hz, A0_5__4Hz};

AmplitudeSineFreqInfo = {'A = 1, f_{sin} = 0.0960Hz', 'A = 1, f_{sin} = 0.1330Hz', 'A = 1, f_{sin} = 0.1664Hz', 'A = 1, f_{sin} = 0.2470Hz', 'A = 1, f_{sin} = 0.5004Hz', 'A = 1, f_{sin} = 1.0164Hz', 'A = 1, f_{sin} = 2.0151Hz', 'A = 1, f_{sin} = 4.0099Hz'};

%%
iter = 0;

plot_iterator = 1;

% two-row table - first row: A (complex num. magnitude), second: angle 
MagnitudePhaseInput = zeros(2,length(sineFrequencies)); %data from input signal 
MagnitudePhaseOutput = zeros(2,length(sineFrequencies)); %from output signal
Gs = zeros(1, length(sineFrequencies)); % gain
fis = zeros(1, length(sineFrequencies)); % phase


ReWe=zeros(1,length(sineFrequencies)); %real part of input signal
ImWe=zeros(1,length(sineFrequencies)); %Imaginary part of output signal

ReWy=zeros(1,length(sineFrequencies));
ImWy=zeros(1,length(sineFrequencies));

for iter=1:length(matrixOfResults);

matrixOfResults{iter}(2,:)=matrixOfResults{iter}(2,:).*1.41;
Nt = size(matrixOfResults{iter}, 2); %number of samples

fp= 81.4;% <- sampling frequency
fn=fp/2; %Nyquist frequency



t=(0:Nt-1)/fp; %time vector needed for DFT

Nf = round(Nt/2+1); %sample corresponding to Nyquist frequency





% for fi=1:length(sineFrequencies)

    freq=sineFrequencies(iter); %calculating DFT for a specific test freq.

    ReWe(iter) = 2/Nt * korL(matrixOfResults{iter}(1,:), cos(2*pi*freq*t)); % real part DFT 
    ImWe(iter) = 2/Nt * korL(matrixOfResults{iter}(1,:), sin(2*pi*freq*t)); %imaginary part DFT

    ReWy(iter) = 2/Nt * korL(matrixOfResults{iter}(2,:), cos(2*pi*freq*t));
    ImWy(iter) = 2/Nt * korL(matrixOfResults{iter}(2,:), sin(2*pi*freq*t));
    
% end



% normalization 
% ReWe(1) = ReWe(1)/2;
% ReWe(Nf) = ReWe(Nf)/2;
% 
% ImWe(1) = ImWe(1)/2;
% ImWe(Nf) = ImWe(Nf)/2;
% 
% ReWy(1) = ReWy(1)/2;
% ReWy(Nf) = ReWy(Nf)/2;
% 
% ImWy(1) = ImWy(1)/2;
% ImWy(Nf) = ImWy(Nf)/2;

MagnitudePhaseInput(1,iter)=sqrt(ReWe(iter)^2 + ImWe(iter)^2); %calculate input sig. Magnitude
MagnitudePhaseInput(2, iter)= rad2deg(atan2(ImWe(iter),ReWe(iter))); %ditto - phase 360 deg / 256 = 1.41

MagnitudePhaseOutput(1,iter)=sqrt(ReWy(iter)^2 + ImWy(iter)^2);
MagnitudePhaseOutput(2,iter)= rad2deg(atan2(ImWy(iter),ReWy(iter)));


end 

Gs = MagnitudePhaseOutput(1,:) ./ MagnitudePhaseInput(1,:);
fis = MagnitudePhaseOutput(2,:) - MagnitudePhaseInput(2,:); 



    figure;
    t = tiledlayout(2,1);
    
    nexttile
    %plot(sineFrequencies,Gs,'b-');
    stem(sineFrequencies,Gs);
    xlim([0 4.5]);
    title("Ch-ka Bode'ego G(f)");
    xlabel('Hz')
    ylabel('G')
    
    nexttile
    %plot(sineFrequencies,fis);
    stem(sineFrequencies,fis);
    xlim([0 4.5]);
    title("Ch-ka Bode'ego dfi(f)");
    xlabel('Hz')
    ylabel('delta fi')
    %title(t, AmplitudeSineFreqInfo);



