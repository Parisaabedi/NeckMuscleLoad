%% plotting Data (Version 2)
% Written by Parisa Abedi, NOV. 2015
% This function plots the data from the simulation, one value can be fixed 
% and the other value will be varied to see the effect
% (you must have the Wd.mat from your simulation)

clear all
close all

% enter 'V' for varying value and the number for fixed value
hr = inputdlg('Head Roll','Please enter Head Roll here',1);
nl = inputdlg('Neck Load','Please enter Neck Load here',1);
s = dir;
if hr{1,1} == 'V'
    NL = str2num(nl{1,1});
    sectionname = strcat('NL=', nl{1,1});
    c = (' Head Roll (deg)');
    xs = [ 30  0 -30];
    xtick = ['30 CCW' '0' '30 CW']; % data and simulations
else
    HR = str2num(hr{1,1});
    sectionname = strcat('HR=', hr{1,1});
    c = (' Neck Load');
    xs = [ +1  0 -1];
    xtick = ['Left' '0' 'Right'];
end

x = [ 0    45    90   135   180   225   270   315]; % target positions (deg)

colour = ['m','c','r','k','g'];
pf1 = zeros(3,1);

s1 = 0;
s2 = 0;

ystdc = zeros(1,3);
ystdci = zeros(3,5);
temp1 = [];

for j = 1:length(s)
    k = strfind(s(j).name, sectionname);
        if k ~= 0      
            st = s(j).name;
            load(st)
            load('Wd.mat') 
            if hr{1,1} == 'V'
              hp = strfind(st, 'HR='); % find the HR position
              np = strfind(st, 'NL=');
              nle = strfind(st, '.mat');
              hr1 = st(hp+3:np-3); % extract the HR value
              cdn = strcat('HR = ', hr1); % current display name
              hr1 = str2num(hr1);
              j1 = find(xs == hr1);
            else
              np = strfind(st, 'NL='); % find the NL position
              nle = strfind(st, '.mat');
              nl = st(np+3:nle-1); % extract the NL value
              cdn = strcat('NL = ', nl); % current display name
              nl = str2num(nl);
              j1 = find(xs == nl);
            end      

            y2p = zeros(8,5);
            y2v = zeros(8,5);
            te = exist('targetp', 'var');
            if ~te
                targetp = target;
            end
            ystdm = ones(size(targetp,1),5*8)*1000;
            for i=1:8
                for i2= 1:5
                    q = targetp(:,i2+5*(i-1))~=1000;
                    y2p(i,i2) = mean(targetp(q,i2+5*(i-1)));   
                    y2v(i,i2) = mean(targetv(q,i2+5*(i-1)));
                    ystdm(q,i2+5*(i-1)) = targetp(q,i2+5*(i-1))- ones(sum(q),1) * y2p(i,i2);      
                end
            end 
            % Reach Error based on Prop. Movement Vector
            figure(1); subplot(1,2,1);   hold on    
            pf1(j1) = plot(x,y2p(:,4) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle','-', 'DisplayName',cdn);    
            plot(x,y2p(:,2) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle',':')
%             plot(x,y2p(:,3) ,'Color',colour(j1),'LineWidth',3.0,...
%                            'Marker','o','MarkerSize',8, 'LineStyle','--')
            figure(1);   subplot(1,2,2);  hold on;   
            plot(x,y2v(:,4) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle','-', 'DisplayName',cdn);    
            plot(x,y2v(:,2) - wd(:,1),'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle',':')
%             plot(x,y2v(:,3) ,'Color',colour(j1),'LineWidth',3.0,...
%                            'Marker','o','MarkerSize',8, 'LineStyle','--')
                       
            figure(2);  subplot(1,2,1); hold on;   
            pf2(j1) = plot(x,y2p(:,5) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle','-', 'DisplayName',cdn);    
            plot(x,y2p(:,1) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle',':')
%             plot(x,y2p(:,3) ,'Color',colour(j1),'LineWidth',3.0,...
%                            'Marker','o','MarkerSize',8, 'LineStyle','--')
            figure(2);  subplot(1,2,2);  hold on;    
            plot(x,y2v(:,5) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle','-', 'DisplayName',cdn);    
            plot(x,y2v(:,1) - wd(:,1) ,'Color',colour(j1),'LineWidth',3.0,...
                           'Marker','o','MarkerSize',8, 'LineStyle',':')
%             plot(x,y2v(:,3) ,'Color',colour(j1),'LineWidth',3.0,...
%                            'Marker','o','MarkerSize',8, 'LineStyle','--')
            % find the SD for each IHP position seperately (ADDED07/2015)
            for i2 = 1:5
                ysdmI = ystdm(:,i2*5+1:(i2+1)*5); 
                ystdmi = ysdmI(:);
                q = ystdmi(:) ~= 1000;
                temp = ystdmi(q);
                ystdci(j1,i2) = std(temp); 
            end  
            
            % find the Standard deviation
            ystdm = ystdm(:);
            q = ystdm(:)~=1000; % without any threshold
%             q = ystdm(:)~=1000 & abs(ystdm(:)) < 22.5; % with threshold
            temp = ystdm(q);
            ystdc(j1) = std(temp); 
            if j ~= 2
                temp1 = [temp1; temp];   
            end 
            
           
        end              
end

figure(2); subplot(1,2,1)  
set(gca,'XTick',x,'FontSize',12,'FontWeight','bold'); 
xlabel('Target Direction(deg)','FontSize',12,'FontWeight','bold')
ylabel('Visual Reach Error(deg)','FontSize',12,'FontWeight','bold')
title(strcat('Effect of changing',c, ' on reach error for Lateral Shift'))
AX = legend(pf1);  
LEG = findobj(AX,'type','text'); set(LEG,'fontsize',14)

figure(2);subplot(1,2,2)
set(gca,'FontSize',12,'FontWeight','bold'); 
xlabel(c,'FontSize',12,'FontWeight','bold')
ylabel('Visual Reach Error(deg)','FontSize',12,'FontWeight','bold')
title(strcat('Effect of changing',c, ' on reach error for Depth Shift'))

figure(1); subplot(1,2,1)  
set(gca,'XTick',x,'FontSize',12,'FontWeight','bold'); 
xlabel('Target Direction(deg)','FontSize',12,'FontWeight','bold')
ylabel('Prop. Reach Error(deg)','FontSize',12,'FontWeight','bold')
title(strcat('Effect of changing',c, ' on reach error for Lateral Shift'))
AX = legend(pf1);  
LEG = findobj(AX,'type','text'); set(LEG,'fontsize',14)

figure(1);subplot(1,2,2)
set(gca,'FontSize',12,'FontWeight','bold'); 
xlabel(c,'FontSize',12,'FontWeight','bold')
ylabel('Prop. Reach Error(deg)','FontSize',12,'FontWeight','bold')
title(strcat('Effect of changing',c, ' on reach error for Depth Shift'))


figure(3)
bar(xs , ystdc(3:-1:1))
set(gca,'FontSize',12,'FontWeight','bold'); 
xlabel(c,'FontSize',12,'FontWeight','bold')
ylabel('Standard Deviation','FontSize',12,'FontWeight','bold')
title(strcat('Effect of changing',c, ' on Standard deviation'))


