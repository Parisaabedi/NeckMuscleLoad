%% Allocentric (PA: March 2017)
% This function performs the required operations for allocentric
% information of our targeted object
% this version includes causality!
function [Amo,Avarmo, Poly] = Allo_coding_gen(OR,OIR, Tid, varv, varu, Sid, shiftA,R,C,IRR)
    % Initialization
%     varv = 0.1; % visual variability
%     varu = 0.2; % uncertainty due to incongruency
%     % Find object positions
%     O = randi(100,10,2); % For the example code!
%     Ro = 6;
%     % Tid = randi(6,1); % Find the target
%     % ot = [O(1:Tid-1,:) ; O(Tid+1:end,:)]; % Use the other points to creat Triangles (Polygons)
%     ot = O(1:Ro,:);
    %% Encoding
%     load('test.mat')
    ot = OR;
    No = length(OR);Noi = length(OIR);
    %%%%%%%%% Step 1 : Find the Barycenter %%%%%%%%%
    % create a polygon (convex or concave)
    v = boundary(ot(:,1),ot(:,2));
%     plot(ot(v,1),ot(v,2)) % Plot the polygon
    Poly = [ot(v,1) ot(v,2)];
    % Find the cetroid
    NR = length(OR);
    NIR = length(OIR);
    mr = R/NR;
    mir = (1-R)/NIR;
    if mir == 0
        mir = 1e-20;
    end
%     [A, bc(1,1),bc(1,2)] = polycenter(Poly(1:end-1,1),Poly(1:end-1,2));
    cm = mr * sum(OR(:,:)) + mir * sum(OIR(:,:));
%     hold on; scatter(cm(1,1),cm(1,2),'r*')
    % CE = polygonCentroid(Poly);
    % hold on; scatter(CE(1,1),CE(1,2),'b+')

    % Transform from CartesianToBaryCenter
%     V = [O' ; ones(1, length(O))];
%     BC = [cm' ; 1];
%     lambda = pinv(V) * BC; % BC = V*lambda
    
    lambda = [ mr*ones(No,1) ; mir*ones(Noi,1)];
    % calculate the variance for BC 
    varX = eye(No + Noi);
    varV = ones(No + Noi,1);
    for i = 1 : No + Noi
        if i < No || i == No
            varX(i,i) = (varv+0.1/mr);
            varV(i) = (varv+0.1/mr);
        else
            varX(i,i) = (varv+0.1/mir);
            varV(i) = (varv+0.1/mir);
        end   
    end
%     varVBC = varX * eye(length(O));
    varBC = lambda' * varX * lambda;

    % Calculate R vectors (The distance of each object from the centroid)
    R = [OR ; OIR] - cm;
    % calculate the variance for R
    varR = varV + varBC;
    %% Decoding
    % Missing object (Target)
    if IRR == 1 % Shift relevant objects
        Nor = [OR(1:Tid-1,:) ; OR(Tid+1:end,:)]; % Other observable objects
        lambdan = [ lambda(1:Tid-1) ; lambda(Tid+1:end)];
        varRn = [varR(1:Tid-1,:) ; varR(Tid+1:end,:)];
        % ******** shift some of the objects *********************
        Ts = Sid ; 
        Nor(Ts,1) = Nor(Ts,1) + C*shiftA;
        Noir = OIR;
    elseif IRR == 0 % shift irrelevant objects
       % ******** shift some of the objects *********************
        Ts = Sid ;
        Noir = OIR;
        Noir(Ts,1) = OIR(Ts,1) + C*shiftA;
        Nor = [OR(1:Tid-1,:) ; OR(Tid+1:end,:)]; % Other observable objects
        lambdan = [ lambda(1:Tid-1) ; lambda(Tid+1:end)];
        varRn = [varR(1:Tid-1,:) ; varR(Tid+1:end,:)];
    elseif IRR == 2 % both relevant and irrelevant objects are shifted!       
        Nor = [OR(1:Tid-1,:) ; OR(Tid+1:end,:)]; % Other observable objects
        lambdan = [ lambda(1:Tid-1) ; lambda(Tid+1:end)];
        varRn = [varR(1:Tid-1,:) ; varR(Tid+1:end,:)];
        Nor = Nor + C(1:length(Nor)) .* shiftA;
        Noir = OIR + C(length(Nor)+1:end) .* shiftA;
    end
    No = [Nor ; Noir ];

    % Find the new BC
    NR = [R(1:Tid-1,:) ; R(Tid+1:end,:)];
    Nbc = No - NR;
    Nbc = round(Nbc,4);
    Nic = zeros(length(Nbc),1); % find the number of congruency and related weights
    for i = 1 : length(Nbc)
        q = Nbc(:,1) == Nbc(i,1) & Nbc(:,2) == Nbc(i,2);
        Nic(i) = sum(q.*lambdan);
    end

    % Calculate the variance for each new BC (observation + adding the r calculation)
    varfnbc = varv + varRn;
    
    % perform causality: chech which ones has more support and run the
    % integration for that one!
    
    [a,b] = max(Nic);
    q = Nbc(:,1) == Nbc(b,1) & Nbc(:,2) == Nbc(b,2);
    % Calculate the final BC by integrating all the Nbc (Bayesian Integration)
    varnBC = (sum(ones(length(varfnbc(q)),1)./varfnbc(q)))^-1;
    nBC = sum(varnBC.*(Nbc(q,:) ./ [varfnbc(q) varfnbc(q)])) ;
%     varnBC = varu * sum(Nic.*lambdan) + varnBC; 

    Amo = nBC + R(Tid,:); % Calculte the position for missing object
    Avarmo = varnBC + varR(Tid); % variance of missing object

    % return the mu,variance for the allo-centric
    % disp(['Variance of missing object is ' num2str(varmO) ' and the position is ' num2str(mO)])
end