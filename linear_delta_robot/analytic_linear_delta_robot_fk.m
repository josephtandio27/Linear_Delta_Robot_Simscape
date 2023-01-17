function [eePos, eeVel, statusFlag] = analytic_linear_delta_robot_fk(actuatorPos, actuatorVel)
    % actuatorPos 1x3 double array [m]
    % actuatorVel 1x3 double array [m/s]
    % eePos 1x3 double array, end effector [m]
    % eeVel 1x3 double array, end effector [m/s]
    % statusFlag = 1 if successful
    R = 0.45;
    r = 0.04;
    a = R*(sqrt(3)/2);
    b = pi/4;
    h = 0.35;
    d = r*(sqrt(3)/2);
    based_plate_thick = R/15;
    bottom_plate_thick = r/15;
    eeL = 0.03;

    v = d-a;
    n11 = actuatorPos(3)^2-actuatorPos(2)^2-2*v*actuatorPos(2)*cos(b)+2*v*actuatorPos(3)*cos(b);
    n12 = 2*sqrt(3)*v+sqrt(3)*actuatorPos(2)*cos(b)+sqrt(3)*actuatorPos(3)*cos(b);
    n13 = actuatorPos(3)*cos(b)-actuatorPos(2)*cos(b);
    n14 = 2*actuatorPos(2)*sin(b)-2*actuatorPos(3)*sin(b);
    
    n21 = actuatorPos(3)^2-actuatorPos(1)^2+2*v*actuatorPos(3)*cos(b)-2*v*actuatorPos(1)*cos(b);
    n22 = sqrt(3)*(v+actuatorPos(3)*cos(b));
    n23 = 2*actuatorPos(1)*cos(b)+actuatorPos(3)*cos(b)+3*v;
    n24 = 2*actuatorPos(1)*sin(b)-2*actuatorPos(3)*sin(b);
    
    n31 = (n11*n22-n12*n21)/(n12*n23-n13*n22);
    n32 = (n12*n24-n14*n22)/(n12*n23-n13*n22);
    
    n41 = (-n13*n31-n11)/n12;
    n42 = (n32*n13-n14)/n12;
    
    A = n32^2+n42^2+1;
    B = 2*n41*n42-2*n31*n32-2*actuatorPos(1)*sin(b)+2*n32*v+2*actuatorPos(1)*n32*cos(b);
    C = actuatorPos(1)^2+v^2-h^2+n31^2+n41^2+2*v*actuatorPos(1)*cos(b)-2*n31*v-2*actuatorPos(1)*n31*cos(b);
    
    Pzp = (-B+sqrt(B^2-4*A*C))/(2*A);
    Pzn = (-B-sqrt(B^2-4*A*C))/(2*A);
    Pyp = n31-Pzp*n32;
    Pyn = n31-Pzn*n32;
    Pxp = n41+Pzp*n42;
    Pxn = n41+Pzn*n42;
    Pp = [Pxp Pyp Pzp+bottom_plate_thick/2+eeL-based_plate_thick/2]';
    Pn = [Pxn Pyn Pzn+bottom_plate_thick/2+eeL-based_plate_thick/2]';

    if ~isreal(Pp) || ~isreal(Pn)
        error('imaginary result')
    elseif Pp(3) >= 0
        statusFlag = 1;
        eePos = Pp;
    elseif Pn(3) >= 0
        statusFlag = 1;
        eePos = Pn;
    else
        error('Cannot find the Forward Kinematic')
    end

    Jy = zeros(3);
    Jy(1,1) = 2*eePos(1);
    Jy(1,2) = 2*eePos(2)-2*v-2*cos(b)*actuatorPos(1);
    Jy(1,3) = 2*eePos(3)-2*sin(b)*actuatorPos(1);
    Jy(2,1) = 2*eePos(1)-sqrt(3)*(v+cos(b)*actuatorPos(2));
    Jy(2,2) = v+2*eePos(2)+cos(b)*actuatorPos(2);
    Jy(2,3) = 2*eePos(3)-2*sin(b)*actuatorPos(2);
    Jy(3,1) = 2*eePos(1)+sqrt(3)*(v+cos(b)*actuatorPos(3));
    Jy(3,2) = v+2*eePos(2)+cos(b)*actuatorPos(3);
    Jy(3,3) = 2*eePos(3)-2*sin(b)*actuatorPos(3);

    Jq = zeros(3);
    Jq(1,1) = -(2*actuatorPos(1)+2*v*cos(b)-2*eePos(2)*cos(b)-2*eePos(3)*sin(b));
    Jq(2,2) = -(2*actuatorPos(2)+2*v*cos(b)+eePos(2)*cos(b)-2*eePos(3)*sin(b)-sqrt(3)*eePos(1)*cos(b));
    Jq(3,3) = -(2*actuatorPos(3)+2*v*cos(b)+eePos(2)*cos(b)-2*eePos(3)*sin(b)+sqrt(3)*eePos(1)*cos(b));
    
    eeVel = Jy\Jq*actuatorVel;
end