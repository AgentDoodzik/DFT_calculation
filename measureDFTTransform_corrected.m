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
MagnitudePhaseInput = cell(length(sineFrequencies),2); %data from input signal 
MagnitudePhaseOutput = cell(length(sineFrequencies),2); %from output signal
Gs = cell(length(sineFrequencies),1); % gain
fis = cell(length(sineFrequencies),1); % phase

for iter=1:length(matrixOfResults);


Nt = size(matrixOfResults{iter}, 2); %number of samples

fp= 81.4;% <- sampling frequency
fn=fp/2; %Nyquist frequency



t=(0:Nt-1)/fp; %time vector needed for DFT

Nf = round(Nt/2+1); %sample corresponding to Nyquist frequency


ReWe=zeros(1,length(sineFrequencies)); %real part of input signal
ImWe=zeros(1,length(sineFrequencies)); %Imaginary part of output signal

ReWy=zeros(1,length(sineFrequencies));
ImWy=zeros(1,length(sineFrequencies));



for fi=1:length(sineFrequencies)

    freq=sineFrequencies(fi); %calculating DFT for a specific test freq.

    ReWe(fi) = 2/Nt * korL(matrixOfResults{iter}(1,:), cos(2*pi*freq*t)); % real part DFT 
    ImWe(fi) = 2/Nt * korL(matrixOfResults{iter}(1,:), sin(2*pi*freq*t)); %imaginary part DFT

    ReWy(fi) = 2/Nt * korL(matrixOfResults{iter}(2,:), cos(2*pi*freq*t));
    ImWy(fi) = 2/Nt * korL(matrixOfResults{iter}(2,:), sin(2*pi*freq*t));
    
end



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

MagnitudePhaseInput{iter,1}=sqrt(ReWe.^2 + ImWe.^2);
MagnitudePhaseInput{iter,2}= atan2(ImWe,ReWe);

MagnitudePhaseOutput{iter,1}=sqrt(ReWy.^2 + ImWy.^2);
MagnitudePhaseOutput{iter,2}= atan2(ImWy,ReWy);

Gs{iter,1} = MagnitudePhaseOutput{iter,1} ./ MagnitudePhaseInput{iter, 1};
fis{iter,1} = MagnitudePhaseOutput{iter,2} - MagnitudePhaseInput{iter,2}; 

% figure(plot_iterator);
% t=tiledlayout(2,1);
% 
% nexttile
% plot(f, MagnitudePhaseInput{iter,1},'r-');
% xlim([0 5]);
% title("Ch-ka Bode'ego syg. we A(f)");
% 
% 
% nexttile
% plot(f, rad2deg(MagnitudePhaseInput{iter,2}),'b-');
% xlim([0 5]);
% title("Ch-ka Bode'ego syg. we fi(f)");
% title(t,AmplitudeSineFreqInfo{iter});
% plot_iterator = plot_iterator+1;
% 
% figure(plot_iterator);
% t=tiledlayout(2,1);
% 
% nexttile
% plot(f, MagnitudePhaseOutput{iter,1},'r-');
% xlim([0 5]);
% title("Ch-ka Bode'ego syg. wy A(f)");
% 
% nexttile
% plot(f, rad2deg(MagnitudePhaseOutput{iter,2}),'b-');
% xlim([0 5]);
% title("Ch-ka Bode'ego syg. wy fi(f)");
% title(t,AmplitudeSineFreqInfo{iter})
% 
% plot_iterator = plot_iterator+1;
% figure(plot_iterator);
% 
% t=tiledlayout(2,1);
% nexttile
% plot(f, Gs{iter,1},'b-');
% xlim([0 5]);
% title("Ch-ka Bode'ego G(f)")
% 
% nexttile
% plot(f,fis{iter,1},'-b');
% xlim([0 5]);
% title("Ch-ka Bode'ego dfi(f)")
% title(t,AmplitudeSineFreqInfo{iter})
% 
% plot_iterator = plot_iterator+1;





end 

%Plotting section (uncomment for actual plotting)
for iter = 1:length(matrixOfResults);
    figure;
    t = tiledlayout(2,1);
    
    nexttile
    %plot(sineFrequencies, Gs{iter, 1},'b-');
    stem(sineFrequencies,Gs{iter,1});
    xlim([0 5]);
    title("Ch-ka Bode'ego G(f)");
    
    nexttile
    %plot(sineFrequencies, fis{iter, 1},'-b');
    stem(sineFrequencies,Gs{iter,1});
    xlim([0 5]);
    title("Ch-ka Bode'ego dfi(f)");
    title(t, AmplitudeSineFreqInfo{iter});
end


