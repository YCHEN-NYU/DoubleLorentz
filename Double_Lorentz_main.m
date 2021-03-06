%% ========================================================================
% FMR Data Fitting main function 
% Author: Yizhang Chen, NYU
% Contact: yc1224@nyu.edu
% Fit FMR data with Double Lorentzian Model

% =========================================================================
% Input data format: *_(frequency)MHz.dat (e. g. para_10000MHz.dat) 
% 1     2    3       4        5      6          7       8        9      10
% Index H(T) S11(Re) S11(Img) S11(A) S11(Phase) S12(Re) S12(Img) S12(A) S12(Phase) 

% 11      12       13     14      
% S22(Re) S22(Img) S22(A) S22(Phase)

% =========================================================================
% Output data format (*.txt)

% frequency(GHz) Hres1 w1 w_lb1 w_ub1 Hres2 w2 w_lb2 w_ub2


% Output figure format (*.png) with the input data file names (e. g. para_1000MHz.png)

%% ========================================================================
% clean up existing data & close all open windows
clear;
close all;

%% ========================================================================
% Give destination folder, outputfile
cd '/Users/yiyi/Desktop/FMR_dataanalysis/f-H/#17/17_P2N1';
outputname='fit_parameters.txt';
% Give indices of field(T), S12_real, S12_imag
indexH = 2;
indexS12real = 7;
indexS12imag = 8;

%% ========================================================================

%%% read files in the current folder
folder = pwd;
% read all data file ended with *MHz.dat in the current folder
files=dir('*MHz.dat');
% sort out the files in natural order
[filenames, index] = sort_nat({files.name});
% open and creat the file for output
outputloc=[folder '/' outputname];
% Open or create new file for writing. Append data to the end of the file.
fidout=fopen(outputloc,'a+');


%% ========================================================================
% starting and ending indices of data files
i=14;
% i_end = length(filenames);
i_end = 15;
while i<=i_end;
%% ========================================================================
% Get frequency(GHz) from file names (read out number between _ and MHz)
str_splitted = strsplit(char(filenames{i}),'_');
strend = str_splitted(end);% Take '*MHz.dat' from data files 'XXX_*MHz.dat'
final_str = strsplit(char(strend),'MHz'); % Get '*' by splitting '*MHz.dat'
frequency = str2double(char(final_str(1)))/1000;% convert MHz to GHz

%% ========================================================================
% load data
fileloc = [folder '/' char(filenames(i))];
fprintf('%s\n',char(filenames(i)));
rawdata = importdata(fileloc);
data = rawdata(2:end-1,:);% get rid of the 1st data point.

%% ========================================================================
% field in Tesla
x = data(:,indexH); 
% Construct complex data y with S12_real and S12_imag
y = complex(data(:,indexS12real),data(:,indexS12imag)); 

[fitpara, fitconfint, S12_plot] = Double_Lorentz(x,y,frequency);

%% ========================================================================
% Type 'y' or 'Y' to save parameters, figures & continue with next data
% files. Type any other characters to re-fit current data.
replyw = input('Want to fit next one?(y/n)', 's');

% if input character is 'y', save fitting parameters & figure
if strcmpi(replyw,'y')    
    % save fitting parameters & confidential intervals (95%)
    Hres1=abs(fitpara(4))*1e4;%%% Change to Oe 
    w1=abs(fitpara(2))*1e4; % in Oe
    w_lb1=min(abs(fitconfint(2,:)))*1e4; % in Oe
    w_ub1=max(abs(fitconfint(2,:)))*1e4; % in Oe
    
    Hres2=abs(fitpara(8))*1e4;%%% Change to Oe 
    w2=abs(fitpara(6))*1e4; % in Oe
    w_lb2=min(abs(fitconfint(6,:)))*1e4; % in Oe
    w_ub2=max(abs(fitconfint(6,:)))*1e4; % in Oe
    
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';% set image size as auto

    saveas(S12_plot,strtok(char(filenames{i}),'.'),'png');

        
    % Save the fitting result 
    % frequency,Hres1,w1,w_lb1,w_ub1,Hres2,w2,w_lb2,w_ub2
    fprintf(fidout,'%10.0f %10.2f %10.3f %10.3f %10.3f %10.2f %10.3f %10.3f %10.3f\n',frequency,Hres1,w1,w_lb1,w_ub1,Hres2,w2,w_lb2,w_ub2);

 % Move on to the next data file
    i=i+1;
end
            close(S12_plot);

end

%% ========================================================================

% Print out ending information & close output file
fprintf('End of fitting\n')
fclose(fidout);


