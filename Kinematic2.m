% This function Calculates the Kinematic and inverse kinematic of the arm
% written by Parisa Abedi, NOV. 2015
function xy = Kinematic2(theta,kin,FI)
if FI == 1 %(Forward Kinematic)
    t1 = theta(1,1,:,:);
    t2 = theta(2,1,:,:);
    l1= kin(1);
    l2 =  kin(2);
    xy(1,1,:,:) = -sin(t1) .* (l1*sin(t2)+l2);
    xy(2,1,:,:) =  cos(t1) .* (l1*sin(t2)+l2);
elseif FI == -1 %(Inverse Kinematic)
    x = theta;    
    l1= kin(1);
    l2 =  kin(2);
    t1 = atan(-x(1,1,:,:)./x(2,1,:,:));
    temp = sqrt(x(1,1,:,:).^2 + x(2,1,:,:).^2) - l2;
    t2 = asin(temp/l1); 
    if ~isreal(t2)
        t2 = real(t2);        
    end
    xy(1,1,:,:) = t1;
    xy(2,1,:,:) = t2;    
else
    error('myApp:argChk', 'Wrong Input')
end

end