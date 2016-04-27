%%%%%%%%%%%%this function basically needs the following data
%Input x,y = y_real+i*y_imag, frequency(GHz)
% Fit double Lorentzian 


function [fitpara, fitconfint, S12_plot] = Double_Lorentz(x,y,frequency)

% ==========================================
% Plot out full S12-H figures

            format long; % set output format as long
            S12_plot = figure();
            set(S12_plot, 'Position', [200, 100, 800, 600])
            set(gcf,'color','w');
            % S12 Real
            subplot(2,1,1);
            plot(x,real(y),'ro','markersize',10);
            ylabel('S12 Real','FontSize',36,'FontWeight','bold') 
            set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
            set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
            set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));
            % S12 Imaginary
            subplot(2,1,2);
            plot(x,imag(y),'bs','markersize',10);
            xlabel('H(T)','FontSize',36,'FontWeight','bold')
            ylabel('S12 Img','FontSize',36,'FontWeight','bold') 
            set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
            set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
            set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));
            % set a common title on the top of the figure
            set(gcf,'NextPlot','add');
            axes; 
            set(gca,'Visible','off'); 
            h = title(['f = ' num2str(frequency) 'GHz'],'fontsize',40,'fontweight','b');
            set(h,'Visible','on');
                        
            %%% set up the region for plotting
            plotxR=0.1;
%             plotyR=0.1;
            
            xlim([min(x)-plotxR*(max(x)-min(x)) max(x)+plotxR*(max(x)-min(x))]);
            
%           ylim([min(real(y))-plotyR*(max(real(y))-min(real(y))) max(y)+plotyR*(max(y)-min(y))]);
            set(gca,'Fontsize',30,'Linewidth',3);
            title(['f=' num2str(frequency) 'GHz'])
           
            %%% To select the region for fitting
            fprintf('Click on [LEFT, RIGHT] Region!\n')
            
            [xleft,~]=ginput(1);
            [xright,~]=ginput(1);
            
            % swap xleft and xright if xleft > xright
            if(xleft > xright)
                xtemp = xright;
                xright = xleft;
                xleft = xtemp;
            end
            
            %%% This loop exclude the data that is outside the region
            %%% selected. 
            ind=zeros(length(x),1);
            for i =1:1:length(x);
                if x(i)<xleft || x(i)>xright
                   ind(i)=i;
                end               
            end
            
            x=x(setdiff(1:length(x),ind));
            y=y(setdiff(1:length(y),ind));
            
            %%% reset the figure and plot with the region just selected. 
            % ==========================================================
            clf(S12_plot,'reset');
            hold on;

            % S12 Real
            subplot(2,1,1);
            h1 = plot(x,real(y),'ro','markersize',10);
            deltax = 0.05*(max(x)- min(x));
            xlim([min(x)-deltax,max(x)+deltax]);
            deltay_real = 0.05*(max(imag(y))- min(imag(y)));
            ylim([min(real(y))-deltay_real,max(real(y))+deltay_real]);
            ylabel('S12 Real','FontSize',36,'FontWeight','bold') 
            set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
            set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
%             set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.1g'));

            % S12 Imaginary
            subplot(2,1,2);
            
            h2 = plot(x,imag(y),'bs','markersize',10);
            xlabel('H(T)','FontSize',36,'FontWeight','bold')
            ylabel('S12 Img','FontSize',36,'FontWeight','bold') 
            set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
            xlim([min(x)-deltax,max(x)+deltax]);
            deltay_imag = 0.05*(max(imag(y))- min(imag(y)));
            ylim([min(imag(y))-deltay_imag,max(imag(y))+deltay_imag]);
            set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
%             set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.1g'));
            
            set(gcf,'NextPlot','add');
            axes; 
            set(gca,'Visible','off'); 
            h = title(['f = ' num2str(frequency) 'GHz'],'fontsize',40,'fontweight','b');
            set(h,'Visible','on');


%%
%%% just checking testx and testy are finite (wich are the plot we want to fit)
            %%% exclude data exclude points from testx,testy through the logic ok_ variable.
            %%% if ok_ is zero we exclude that point
            
            ok_ = isfinite(x) & isfinite(real(y)) & isfinite(imag(y)); %&ok1&ok2&ok3; 
            x=x(ok_); 
            y=y(ok_);
            
% ==================================
% store data to be fitted inside xdata (n,1) and ydata (n,2)
            xdata = x;
            ydata = [real(y),imag(y)];
% ==================================
       
            
            %%
%Initialization points construction
% theta = n*pi, symmetric; theta = (n+1/2)*pi, asymmetric

% ==================================
% 1st Peak LEFT 
            fprintf('Click on the LEFT of 1st Peak!\n');
            [x_left_1st,~] = ginput(1);
% find the closest point by x_1st position
            dtemp = abs(xdata-x_left_1st*ones(size(xdata)));
            [~,i_left_1st] = min(dtemp);
            y_left_1st =ydata(i_left_1st,:);
            clear dtemp;
  % ==================================
% 1st Peak Right 
            fprintf('Click on the RIGHT of 1st Peak!\n');
            [x_right_1st,~] = ginput(1);
% find the closest point by x_1st position
            dtemp = abs(xdata-x_right_1st*ones(size(xdata)));
            [~,i_right_1st] = min(dtemp);
            y_right_1st =ydata(i_right_1st,:);
            clear temp;
   % ==================================
% 1st Peak Center 
            fprintf('Click on the CENTER of 1st Peak!\n');
            [x_1st,~] = ginput(1);
% find the closest point by x_1st position
            dtemp = abs(xdata-x_1st*ones(size(xdata)));
            [~,i_1st] = min(dtemp);
            y_1st =ydata(i_1st,:);     
            clear dtemp;
                   
% ==================================
% 2nd Peak LEFT 
            fprintf('Click on the LEFT of 2nd Peak!\n');
            [x_left_2nd,~] = ginput(1);
% find the closest point by x_1st position
            dtemp = abs(xdata-x_left_2nd*ones(size(xdata)));
            [~,i_left_2nd] = min(dtemp);
            y_left_2nd =ydata(i_left_2nd,:);
            clear dtemp;
  % ==================================
% 2nd Peak Right 
            fprintf('Click on the RIGHT of 2nd Peak!\n');
            [x_right_2nd,~] = ginput(1);
% find the closest point by x_1st position
            dtemp = abs(xdata-x_right_2nd*ones(size(xdata)));
            [~,i_right_2nd] = min(dtemp);
            y_right_2nd =ydata(i_right_2nd,:);
            clear temp;
   % ==================================
% 2nd Peak Center 
            fprintf('Click on the CENTER of 2nd Peak!\n');
            [x_2nd,~] = ginput(1);
% find the closest point by x_1st position
            dtemp = abs(xdata-x_2nd*ones(size(xdata)));
            [~,i_2nd] = min(dtemp);
            y_2nd =ydata(i_2nd,:);     
            clear dtemp;
                   
%==================================
% %%% Construct starting parameters (11 of them)
% Real 
              % 1st peak
              st_a11 = abs(0.5*(2*y_1st(:,1) - y_left_1st(:,1) - y_right_1st(:,1)));
              st_w11 = x_right_1st-x_left_1st;
              st_theta11 = 0;% start with real symmetric and imaginary asymmetric
              st_x11 = x_1st;

              % 2nd peak
              st_a12 = abs(0.5*(2*y_2nd(:,1) - y_left_2nd(:,1) - y_right_2nd(:,1)));
              st_w12 = x_right_2nd-x_left_2nd;
              st_theta12 = 0;% start with real symmetric and imaginary asymmetric
              st_x12 = x_2nd;
% Imaginary 
              % 1st peak
              st_a21 = abs(0.5*(2*y_1st(:,2) - y_left_1st(:,2) - y_right_1st(:,2)));
%               st_w21 = st_w11;
%               st_theta21 = 0;% start with real symmetric and imaginary asymmetric
%               st_x21 = x_1st;

              % 2nd peak
              st_a22 = abs(0.5*(2*y_2nd(:,2) - y_left_2nd(:,2) - y_right_2nd(:,2)));
%               st_w22 = st_w21;
%               st_theta22 = 0;% start with real symmetric and imaginary asymmetric
%               st_x22 = x_2nd;
              
              
              st_c10 = mean(ydata(:,1));
              st_c11= (ydata(length(ydata),1)-ydata(1,1))/(xdata(length(xdata))-xdata(1));
              st_c12 = 0;
              
              st_c20 = mean(ydata(:,2));
              st_c21= (ydata(length(ydata),2)-ydata(1,2))/(xdata(length(xdata))-xdata(1));
              st_c22 = 0;

              st_ = [st_a11, st_w11, st_theta11,st_x11, st_a12,st_w12, st_theta12,st_x12, st_c10,st_c11,st_c12,st_c20,st_c21,st_c22,st_a21,st_a22];
%==================================
% set fitting options
options = optimset('Tolx',1e-10,'TolFun',1e-10,'DiffMinChange', 1e-16,'MaxIter',15000,'MaxFunEvals',2000000);

% define lower bounds and upper bounds of fitting parameters
x0_lb = [-1,0,-2*pi,0,-1,0,-2*pi,0,-1,-0.5,-0.1,-1,-0.5,-0.1,-0.1,-0.1];
x0_ub = [1,0.5,2*pi,1,1,0.5,2*pi,1,1,0.5,0.1,1,0.5,0.1,0.1,0.1];
%%
%==================================
% fit curve by least square non-linear model defined in cplx_fun.m
% p: fitting parameters
[p,resnorm,residuals,exitflag,output,lambda,jacobian] = lsqcurvefit(@Double_Lorentz_fun,st_,xdata,ydata,x0_lb,x0_ub,options);


% ==================================
% Plot out fitting curves for both real and imaginary parts
xmesh = linspace(min(x),max(x),1000);
yfit1 = p(1)*(p(2)*cos(p(3))+(xmesh-p(4))*sin(p(3)))./(p(2)^2+(xmesh-p(4)).^2)+...
    + p(5)*(p(6)*cos(p(7))+(xmesh-p(8))*sin(p(7)))./(p(6)^2+(xmesh-p(8)).^2) + ...
    p(9)+p(10)*xmesh+p(11)*xmesh.^2;

yfit2  = p(15)*(p(2)*cos(p(3)+pi/2)+(xmesh-p(4))*sin(p(3)+pi/2))./(p(2)^2+(xmesh-p(4)).^2)+...
         + p(16)*(p(6)*cos(p(7)+pi/2)+(xmesh-p(8))*sin(p(7)+pi/2))./(p(6)^2+(xmesh-p(8)).^2) + ...
         p(12)+p(13)*xmesh+p(14)*xmesh.^2;

subplot(2,1,1);
h3 = line(xmesh,yfit1,'linewidth',2,'color','r');
lgd1 = legend([h1,h3],'real','real-fit','location','northeast');
set(lgd1,'fontsize',14);
% set(gca, 'YTick', Y, 'YTickLabel', sprintf('%0.3f|', y));
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
% set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));

subplot(2,1,2);
h4 = line(xmesh,yfit2,'linewidth',2,'color','b');
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
% set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));

lgd2 = legend([h2,h4],'imag','imag-fit','location','northeast');
set(lgd2,'FontSize',14);

% set a common title for the figure;
set(gcf,'NextPlot','add');
axes; 
set(gca,'Visible','off'); 
h5 = title(['f = ' num2str(frequency) 'GHz'],'fontsize',40,'fontweight','b');
set(h5,'Visible','on');

% ==================================
% return fitting parameters and 95% confindent intervals
fitpara = p;
fitconfint = nlparci(p,residuals,'jacobian',jacobian);

end

