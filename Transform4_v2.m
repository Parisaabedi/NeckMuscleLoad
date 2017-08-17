% This function performs the Spatial Transformation (PA: NOV.2015)
% This version performs all the multiplication using mulitdimensional
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input
% xyv : Visual Information (Cartesian,Gaze-Centered)
% thp : Proprioceptive information (Joint Angles, Shoulder-Centered)
% E : Eye in Head
% H : Head in Body
% Sh : Shoulder in Body
% Bs : Body in Space
% kin : Kinematics of arm
% CN : Constant transformation noise

% output
% VC: visual and prop. information in Visual Coordinate
% PC: visual and prop. information in Prop. Coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [thv,xyp] = Transform4_v2(xyv,thp,E,H,Sh,Bs,EHL,Depth,kin,CN)

% consider the linkage in the translation
depth = Depth - EHL;

%%Shoulder-to-Retina
p = size(xyv);
E = ones(1,1,8*5,p(4)) * E;
if size(H) == 1
    H = ones(1,1,8*5,p(4)) * H;
end
% Gaze coordinate 
% Forward Kinematic
xyp = Kinematic2(thp,kin,1); % Prop. information in Cartesian (shoulder-Centered)

% Shoulder Frame, based on Gravity
if Bs
    T = [cosd(Bs) sind(Bs); -sind(Bs) cosd(Bs)];
    for i = 1:size(xyp,1)
        xyp(i,:) = T * xyp(i,:)';
    end
end

% Neck-Centered (transform it to the middle of the body)
xyp(1,1,:,:) = xyp(1,1,:,:) + Sh;

% Head Centered
xyp(2,1,:,:) = xyp(2,1,:,:) - depth;

% Head-Centered (Compensate for head rotation) (CW of RF means CCW of Vector)
T = zeros(2,2,5*8,p(4));
T(1,:,:,:) = [cosd(-H) sind(-H)];
T(2,:,:,:) = [-sind(-H) cosd(-H)]; 
xyp = mmx('Mult',T,xyp);

% Eye-Head- Translation
xyp(2,1,:,:) = xyp(2,1,:,:) - EHL;  

% Retinal-Centered (eye in head rotation)
T = zeros(2,2,5*8,p(4));
T(1,:,:,:) = [cosd(-E) sind(-E)];
T(2,:,:,:) = [-sind(-E) cosd(-E)]; 
xyp = mmx('Mult',T,xyp);

% Constant Transformation Noise
xyp = xyp + CN;

%% Retina-to-Shoulder

% Eye-Centered (eye-head related rotation)
T = zeros(2,2,5*8,p(4));
T(1,:,:,:) = [cosd(E) sind(E)];
T(2,:,:,:) = [-sind(E) cosd(E)]; 
xyv = mmx('Mult',T,xyv);

% Eye-Head- Translation
xyv(2,1,:,:) = xyv(2,1,:,:) + EHL;

% Head-Centered
T = zeros(2,2,5*8,p(4));
T(1,:,:,:) = [cosd(H) sind(H)];
T(2,:,:,:) = [-sind(H) cosd(H)]; 
xyv = mmx('Mult',T,xyv);

% Gaze-Centered(Bring it to the body centered, neck!!)
xyv(2,1,:,:) = xyv(2,1,:,:) + depth;

% Shoulder-Centered
xyv(1,1,:,:) = xyv(1,1,:,:) - Sh;

% Shoulder Frame, based on Gravity
if Bs
    T = [cosd(-Bs) -sind(-Bs); sind(-Bs) cosd(-Bs)];
    for i = 1:size(xyv,1)
        xyv(i,:) = T * xyv(i,:)';
    end
end

% Constant Transformation Noise
xyv = xyv + CN;

% Invers Kinematic
thv = Kinematic2(xyv,kin,-1); 


end