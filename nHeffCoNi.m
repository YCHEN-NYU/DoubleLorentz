clear;
close all;
colortable = ['r','b','c','k','g','m','r','b','c','k','g','m','r','b','c','k','g','m'];
markertable = ['o','s','v','^','o','s','v','^','o','s','v','^','o','s','v','^'];

nHeff =[
4	-0.124
5	-0.1665
6	-0.205
4	-0.1115
5	-0.1685
6	-0.208
4	-0.1095
5	-0.1465
6	-0.1835
4	-0.09715
5	-0.118
6	-0.207
4	-0.08445
5	-0.139
6	-0.196
4	-0.07845
5	-0.1485
6	-0.187
];

p = zeros(1,6);
figure();

for i =1:6
    for j = 1:3
    nH = nHeff((i-1)*3+1:(i-1)*3+j,:);
    end
    p(i) = plot(nH(:,1),nH(:,2),'color',colortable(i),'marker',markertable(i),'markersize',10);
    set(gca,'fontsize',36);

    hold on;
    
end


legend('15','16','17','18','19','20');