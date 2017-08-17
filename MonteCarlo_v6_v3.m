% Monte Carlo Simulation 
% Written by Parisa Abedi, NOV 2014
% It calculated reach movement planning to different visual targets during
% different experimental conditions 
% (Head roll = +30,0,-30, Neck load= right, no load, left, 
%  and IHP = up, down, center, left, and right)

%% Parameters
% av = weight of visual information in visual coordinate
% ap = weight of visual inforamtion in proprioceptive coordinate
% xv = visual position 
% xp = actual proprioceptive position 

%% Initialization
% close all
clear all
sn = 1000; % Simulation Number
l1 = 300; % upper arm length in milimeter
l2 = 450; % lower arm length in milimeter
kin = [l1,l2];
tts = [-250 300]; % Visual Target to Shoulder Transformation Function
ttr = 300; % Visual Target to Retina Transformation
EHL = 130*cosd(71); % Eye-Head linkage 
beta = 1.08; % Head angle estimation gain factor
%(participants rolled their head forward and rested int in the shoulder rest which resulted in 71deg)

% different initial hand positions [yx0xy]
ihp(:,1,1) = [0 -25]'; 
ihp(:,1,2) = [-25 0]';
ihp(:,1,3) = [0 0]';
ihp(:,1,4) = [+25 0]';
ihp(:,1,5) = [0 +25]';
hp = repmat(ihp, [1,1,8,sn]); % create the 3D matrix to make efficient matrix calculation

% Visual Targets distributed on a circle evenly
x = [ 0    45    90   135   180   225   270   315];
tarxp = zeros(1,1,40);
for i = 1:8
    tarxp(1,1,(i-1)*5+1:i*5) = x(i);
end
tarx = repmat(tarxp,[1,1,1,sn]);
R = 100;
xyt(1,1,:,:)= R*cosd(tarx);
xyt(2,1,:,:)= R*sind(tarx);
Target = xyt;
Targeti = Target;
clear x; clear xyt; 

% Visual and Proprioceptive hand postion
xv = zeros(2,1,5*8,sn);
xp = hp;
xip = xp; 

%% Noise of Different Modalities (Direven from fitting the model to data)
Pvar = (60*6.44e-6) ; % Proprioceptive noise
Vvar = (55*0.30) ;    % Visual noise
Tfvar = (90*0.597) ;  % fixed transformation noise
Tvar =  (60*2.46e-3) ;% variable transfroamtion noise
sigmah0 = (60*0.08);  % noise in straight head sensation
MPvar = (27e0); % finalmovement variability induced due to motor movements
k = 25; % contribution of visual&bestibular info. vs. neck muscle info. (sigmavv/sigmanm)

%% Adding Noise to the System
% Visual Noise 
xvv = xv + random('norm',0,sqrt(Vvar),[2,1,5*8,sn]); 
% Proprioceptive noise
% Calculate the joint angles (proprioceptive inforamiton are coded in joint angle coordinates)
CN = zeros(2,1,5*8,sn);
[tp,xyp] = Transform4_v2(xp,CN,0,0,-tts(1),0,EHL,tts(2),[l1,l2],CN); 
% Adding Noise
thp(1,1,:,:) = tp(1,1,:,:) + random('norm',0,sqrt(Pvar),[1,1,40,sn]);
thp(2,1,:,:) = tp(2,1,:,:) + random('norm',0,sqrt(Pvar),[1,1,40,sn]);
% % (tetha -> x) 
[thv,xpp] = Transform4_v2(xvv,thp,0,0,-tts(1),0,EHL,tts(2),[l1,l2],CN);
 % Constant Transformaiton Noise 
CN = random('norm',0,sqrt(Tfvar),[2,1,5*8,sn]);

%% Estimating Head Angle
nh = Tvar; % Head-Roll-Dependent Noise
nl = -1:1:1; % Different Neck Load
ha = 30 * nl; % Different Head Angle

sigmah = nh*abs(ha)+ sigmah0; % associated noise with head angle estimation
uvv = ha; up = ha;

sigmavv = sigmah * (1 +k)/k; % asscociated noise with visual/vestibular informaiton 
sigmap = (1./sigmah - 1./sigmavv).^-1; % asscociated noise with neck muscle informaiton

%% Different Conditions
for ij = 1:length(nl)
    for ic = 1:length(ha) 
            
        % load Control Value
        if ij == 2 && ic == 2
            control = 1;  
        else
            control = 0;
        end
        if control == 0
            load('Wd.mat')
        end
        
        % Different Neck Load and Head Angle
        hr = num2str(ha(ic));
        nlc = num2str(nl(ij));
        sectionname = strcat('HR= ', hr, ' & ',' NL= ', nlc);
        
        % Head angle Estimation
        NL = 1;
        if ic == 2
            HR = 0;
            bi = 1;
        else
            HR = 1;
            bi = 0;
        end
        upt = HR * up(ic) + NL * up(ij); % estimated neck muslce heand angle sense
        sigmapt = mod(abs(upt),29) * sigmap(1) + (upt == 0) * sigmap(2);
        sigmaht = (1/sigmapt + 1/sigmavv(ic))^-1;        
        hvv = random('norm',uvv(ic),sqrt(sigmavv(ic)),[1,1,5*8,sn]);
        hp = random('norm',upt,sqrt(sigmapt),[1,1,5*8,sn]);
        avv = sigmaht / sigmavv(ic); % visual/vestibular weight for heand angle estimation
        eha = avv * hvv + (1-avv) * hp; % MSI for head angle estimation
        eh = reshape(eha,[],1);        
        th = beta * eha; 
        hae = reshape(th,[],1);
        HAE = mean(hae);
        HAV = var(hae-HAE);      
     
        % Transformations
        % point Transformation for the IHPs
        % Noisy head angle Estimation
        xyv = xvv;  
        E = -ha(ic); Bs = 0; Sh = -tts(1); depth = tts(2); kin = [l1 l2];
        [thvn,xypn] = Transform4_v2(xyv,thp,E,th,Sh,Bs,EHL,depth,kin,CN);  
              
                        
        %% Calculating Movement Vector (Visual Coordinate)
        
        % Multi-sensory Integration
            % Calculate the weight for integration  
        
        covv = cell(5,1); 
        covp = cell(5,1);
        av = cell(5,1);      
        for i2 = 1 : 5  
            q = xp(1,1,:,:) == ihp(1,1,i2) & xp(2,1,:,:) == ihp(2,1,i2);             
            tempv = reshape(xvv(:,:,q),2,[]);
            tempp = reshape(xypn(:,:,q),2,[]);
            covvt = cov(tempv'); 
            covpt = cov(tempp'); 
            covtot = inv(covvt\eye(2,2) + covpt\eye(2,2));
            av{i2,1} = 0.5 * (covvt\covtot - covpt\covtot + eye(2,2)); 
            covv{i2,1} = covvt;
            covp{i2,1} = covpt;
        end
        
        Ehpv = zeros(2,1,40,sn);
        reverse = repmat(eye(2),[1,1,40,sn]);
        for i2 = 1 : 5
            q = xp(1,1,:,:) == ihp(1,1,i2) & xp(2,1,:,:) == ihp(2,1,i2);
            avtemp = repmat(av{i2,1},[1,1,40,sn]);
            Ehpv(:,:,q) = mmx('MULT',avtemp(:,:,q),xyv(:,:,q)) + mmx('MULT',(reverse(:,:,q)-avtemp(:,:,q)),xypn(:,:,q));       
        end    
        MV = Target - Ehpv;
        clear avtemp; 

        %% Calculating Initial Joint Angle         
        % Multi-sensory Integration (estimating Initial Joint Angles)
        copv = cell(5,1);
        copp = cell(5,1);
        ap = cell(5,1); 
        for i2 = 1:5
            q = xp(1,1,:,:) == ihp(1,1,i2) & xp(2,1,:,:) == ihp(2,1,i2);             
            tempv = reshape(thvn(:,:,q),2,[]);
            tempp = reshape(thp(:,:,q),2,[]);
            copvt = cov(tempv'); 
            coppt = cov(tempp');             
            covtot = inv(copvt\eye(2,2) + coppt\eye(2,2));
            ap{i2,1} = 0.5 * (copvt\covtot - coppt\covtot + eye(2,2)); 
            copv{i2,1} = copvt;
            copp{i2,1} = coppt;
        end             
         
        ti = zeros(2,1,40,sn);
        for i2 = 1 : 5
            q = xp(1,1,:,:) == ihp(1,1,i2) & xp(2,1,:,:) == ihp(2,1,i2);
            aptemp = repmat(ap{i2,1},[1,1,40,sn]);
            ti(:,:,q) = mmx('MULT',aptemp(:,:,q),thvn(:,:,q)) + mmx('MULT',(reverse(:,:,q)-aptemp(:,:,q)),thp(:,:,q));               
        end
        clear aptemp; clear reverse;
     

        % Inverse Kinematic
        % Velocity Destortion Matrix (VDM)
        Jt = zeros(2,2,5*8,sn);Jti = zeros(2,2,5*8,sn);
        Jt(1,:,:,:) = [-cos(tp(1,1,:,:)).*(l1*sin(tp(2,1,:,:))+l2)  -l1*sin(tp(1,1,:,:)).*cos(tp(2,1,:,:))]; 
        Jt(2,:,:,:) = [-sin(tp(1,1,:,:)).*(l1*sin(tp(2,1,:,:))+l2)   l1*cos(tp(1,1,:,:)).*cos(tp(2,1,:,:))];
        Jti(1,:,:,:) = [-cos(ti(1,1,:,:))./(l1*sin(ti(2,1,:,:))+l2) -sin(ti(1,1,:,:))./(l1*sin(ti(2,1,:,:))+l2)]; 
        Jti(2,:,:,:)  = [-sin(ti(1,1,:,:))./(l1*cos(ti(2,1,:,:)))   cos(ti(1,1,:,:))./(l1*cos(ti(2,1,:,:)))];
        
        DMV = mmx('MULT',Jt,Jti);
               
        % MV vector transformation
        eHa = ha(ic) - th;
        T = zeros(2,2,5*8,sn);
        T(1,:,:,:) = [cosd(-eHa) sind(-eHa)];
        T(2,:,:,:) = [-sind(-eHa) cosd(-eHa)]; 
        MVp = mmx('Mult',T,MV);
        
        % Calculate the Executing Movement angle x' = J(tetha).J(tethai). MV
        % In Proprioceptive Coordinate
        tethapp = zeros(1,1,5*8,sn);
        x = mmx('Mult',DMV,MVp);
        
        xc = x;
        x = reshape(x,2,[]);
        x = x';
        tethap = (180/pi) * atan2(x(:,2),x(:,1));
        tethap = reshape(tethap', sn,[]);
        tethap = reshape(tethap,[],1);
        tethavar = random('norm',0,sqrt(MPvar),size(tethap));
        tethap = tethap + tethavar;
        
        r = sqrt(xc(:,1).^2+xc(:,2).^2);

        % Calculate Disired Movement Angle    
        % Visual Movement
        Targetj = reshape(Targeti,2,[]);
        Targetj = Targetj'; 
        ypxv = (180/pi) * atan2(Targetj(:,2),Targetj(:,1));
        dav = ypxv(:);

        % Proprioceptive Movement
        xipj = reshape(xip,2,[]);  
        xipj = xipj';
        Targetp = Targetj - xipj;        
        ypxp = (180/pi) * atan2(Targetp(:,2),Targetp(:,1));
        dap = ypxp(:);

        % Calculate Error
        Errorv = mod(dav,360) - mod(tethap,360);
        k = abs(Errorv) > 100;
        Errorv(k) = dav(k) - tethap(k);
        
        k = abs(Errorv) < 1e-4;
        Errorv(k) = 0;
        predictedyv = Errorv;

        Errorp = mod(dap,360) - mod(tethap,360);
        k = abs(Errorp) > 100;
        Errorp(k) = dap(k) - tethap(k);
        
        k = abs(Errorp) < 1e-4;
        Errorp(k) = 0;
        predictedyp = Errorp;

        targetv = ones(sn,5*8)*1000;
        targetp = ones(sn,5*8)*1000;
        p = (1:sn)';
        for i = 1:8
            for i2 = 1:5
                targetv(:,(i-1)*5+i2) = predictedyv(((i-1)*5+i2)+(p-1)*40);
                targetp(:,(i-1)*5+i2) = predictedyp(((i-1)*5+i2)+(p-1)*40);
            end
        end
        if control == 1
            y2 = zeros(8,5);
            y2p = zeros(8,5);
            for i=1:8
                for i2= 1:5
                    q = targetv(:,i2+5*(i-1))~=1000;
                    y2(i,i2) = mean(targetv(q,i2+5*(i-1))); 
                    y2p(i,i2) = mean(targetp(q,i2+5*(i-1)));
                end
            end
            wd(:,1) = y2(:,3);
            save('Wd.mat','wd')
        end
        
        save([sectionname '.mat'],'targetv','targetp','av','ap','covv','covp','copv','copp','HAE','HAV')
    end
end
