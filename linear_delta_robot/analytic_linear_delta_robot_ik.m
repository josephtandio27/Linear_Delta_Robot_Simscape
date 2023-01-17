function [actuatorPosDes, statusFlag] = analytic_linear_delta_robot_ik(eePosDes)
    % eePosDes 1x3 double array [m], end effector position destination
    % actuatorPosDes 1x3 double array, actuator position destination [m]
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

    eePosDes(3) = eePosDes(3)-bottom_plate_thick/2-eeL+based_plate_thick/2;
    
    Q1 = eePosDes(2)-d+a;
    G11 = Q1^2+eePosDes(1)^2+eePosDes(3)^2-h^2;
    G12 = -2*eePosDes(3)*sin(b)-2*Q1*cos(b);
    l1p = (-G12+sqrt(G12^2-4*G11))/2;
    l1n = (-G12-sqrt(G12^2-4*G11))/2;
    
    Q21 = eePosDes(1)+(sqrt(3)/2)*(a-d);
    Q22 = eePosDes(2)+(d-a)/2;
    G21 = Q21^2+Q22^2+eePosDes(3)^2-h^2;
    G22 = -sqrt(3)*Q21*cos(b)+Q22*cos(b)-2*eePosDes(3)*sin(b);
    l2p = (-G22+sqrt(G22^2-4*G21))/2;
    l2n = (-G22-sqrt(G22^2-4*G21))/2;
    
    Q31 = eePosDes(1)+(sqrt(3)/2)*(d-a);
    Q32 = eePosDes(2)+(d-a)/2;
    G31 = Q31^2+Q32^2+eePosDes(3)^2-h^2;
    G32 = sqrt(3)*Q31*cos(b)+Q32*cos(b)-2*eePosDes(3)*sin(b);
    l3p = (-G32+sqrt(G32^2-4*G31))/2;
    l3n = (-G32-sqrt(G32^2-4*G31))/2;

    lp = [l1p l2p l3p]';
    ln = [l1n l2n l3n]';

    if ~isreal(lp) || ~isreal(ln)
        error('imaginary result')
    elseif all(ln>=0) && all(ln<=3.5)
        statusFlag = 1;
        actuatorPosDes = ln;
    elseif all(lp>=0) && all(lp<=3.5)
        statusFlag = 1;
        actuatorPosDes = lp;
    else
        error('Cannot find the Inverse Kinematic')
    end
end