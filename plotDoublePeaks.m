clear;
% ==================
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];
% ==================

folder = pwd;
files=dir('*2nd.txt');
[filenames, index] = sort_nat({files.name});% sort out the files in natural order

fff=textscan(filenames{1},'%[^_] %[_] %[^_] %[_] %[^_] %[_] %[^MHz]');
scanname=[cell2mat(fff{1}(1)) cell2mat(fff{2}) cell2mat(fff{3}) cell2mat(fff{4}) cell2mat(fff{5}) cell2mat(fff{6})];
outputname='Heff_2nd.txt';

%%% open and creat the file for output
outputloc=[folder '/' outputname];
fidout=fopen(outputloc,'a+');

i=1;
len_files = length(filenames);
thickness  =[2.67, 2.33, 2.06,1.86,1.69];% nm 

while i<=len_files;
data=importdata(filenames{i});

% frequency,Hres1,w1,w_lb1,w_ub1,Hres2,w2,w_lb2,w_ub2
f = data(:,1);
H1 = data(:,2)/10000;
w1 = data(:,3);
H2 = data(:,6)/10000;
w2 = data(:,8);

fig = figure();
figure(fig);
axes1 = axes('Parent',fig,'FontSize',32);

set(fig, 'Position', [200, 100, 800, 600])

Hmesh = linspace(0,1,1000);

fitP1 = polyfit(H1,f,1);
ffit1 = fitP1(1)*Hmesh+fitP1(2);
Heff1 = fitP1(2)/fitP1(1);

fitP2 = polyfit(H2,f,1);
ffit2 = fitP2(1)*Hmesh+fitP2(2);
Heff2 = fitP2(2)/fitP2(1);

% Save Heff1 & Heff2
fprintf(fidout,'%2.3f %2.3f %2.3f\n',thickness(i),Heff1,Heff2);

h1 = plot(H1,f,'ro','MarkerSize',20);
hold on;
h2 = line(Hmesh,ffit1,'linewidth',5,'color','r');

h3 = plot(H2,f,'bs','MarkerSize',20);

h4 = line(Hmesh,ffit2,'linewidth',5,'color','b');

ylabel('f(GHz)','FontSize',32);
xlabel('H_{res}(T)','FontSize',32)
ylim([0,25]);
xlim([0,1]);
title('Hres - f','fontsize',36);
set(gca,'fontsize',32);
saveas(fig,strtok(char(filenames{i}),'.'),'png')
close(fig);
i = i+1;% move on the the next index

end
