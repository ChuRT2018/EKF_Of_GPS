T = DataShipStruct(2).GPSTime - DataShipStruct(1).GPSTime;

Vx = (ShipCoordinate(2,1) - ShipCoordinate(1,1))/T;
Vy = (ShipCoordinate(2,2) - ShipCoordinate(1,2))/T;
Vz = (ShipCoordinate(2,3) - ShipCoordinate(1,3))/T;
Vt = (ShipCoordinate(2,4) - ShipCoordinate(1,4))/T;

Ax = ((ShipCoordinate(3,1) - ShipCoordinate(2,1))/T - (ShipCoordinate(2,1) - ShipCoordinate(1,1))/T) / T;
Ay = ((ShipCoordinate(3,2) - ShipCoordinate(2,2))/T - (ShipCoordinate(2,2) - ShipCoordinate(1,2))/T) / T;
Az = ((ShipCoordinate(3,3) - ShipCoordinate(2,3))/T - (ShipCoordinate(2,3) - ShipCoordinate(1,3))/T) / T;
At = ((ShipCoordinate(3,4) - ShipCoordinate(2,4))/T - (ShipCoordinate(2,4) - ShipCoordinate(1,4))/T) / T;

X_EV_J(1,1) = ShipCoordinate(1,1) + Vx * T + 0.5*Ax*T*T;
X_EV_J(4,1) = ShipCoordinate(1,2) + Vy * T + 0.5*Ay*T*T;
X_EV_J(7,1) = ShipCoordinate(1,3) + Vz * T + 0.5*Az*T*T;
X_EV_J(10,1) = ShipCoordinate(1,4) + Vt * T + 0.5*At*T*T;

X_EV_J(2,1) = Vx + Ax*T;
X_EV_J(5,1) = Vy + Ay*T;
X_EV_J(8,1) = Vz + Az*T;
X_EV_J(11,1) = Vt + At*T;

X_EV_J(3,1)  = Ax;
X_EV_J(6,1)  = Ay;
X_EV_J(9,1)  = Az;
X_EV_J(12,1) = At;

Dx_ev_J = eye(12) * 100;

De=diag([0.01,0.01,0.01,0.01]);

Recursive_Updating_Result = zeros(12,length(ShipCoordinate));
Recursive_Updating_Result(:,1) = X_EV_J;
%DX(:,1:12) = eye(12) * 4;
vList = zeros(length(ShipCoordinate),1);
T = 0.5;
    diagonal = [
        1 T 0.5*T*T
        0 1 T
        0 0 1
        ];
    Psi = blkdiag(diagonal ,diagonal ,diagonal ,diagonal);
    
   Tell = [
        0.5*T*T 0 0 0 
        T 0 0 0
        0 0 0 0
        0 0.5*T*T 0 0
        0 T 0 0
        0 0 0 0
        0 0 0.5*T*T 0
        0 0 T 0
        0 0 0 0
        0 0 0 0.5*T*T
        0 0 0 T
        0 0 0 0
        ];

for i = 2:(length(ShipCoordinate)+1)
    COUNT = DataShipStruct(i-1).Num;
%     if i <= length(ShipCoordinate)
%         T = DataShipStruct(i).GPSTime - DataShipStruct(i-1).GPSTime;
%     end
%     diagonal = [
%         1 T 0.5*T*T
%         0 1 T
%         0 0 1
%         ];
%     Psi = blkdiag(diagonal ,diagonal ,diagonal ,diagonal);
%     
%     Tell = [
%         T 0.5*T*T   0 0         0 0         0 0
%         1 T         0 0         0 0         0 0
%         0 1         0 0         0 0         0 0
%         0 0         T 0.5*T*T   0 0         0 0
%         0 0         1 T         0 0         0 0
%         0 0         0 1         0 0         0 0
%         0 0         0 0         T 0.5*T*T   0 0
%         0 0         0 0         1 T         0 0
%         0 0         0 0         0 1         0 0
%         0 0         0 0         0 0         T 0.5*T*T
%         0 0         0 0         0 0         1 T
%         0 0         0 0         0 0         0 1
%         ];
    
    for l = 1:COUNT
        X_EV_J0 =  Psi * X_EV_J;
        Dx_ev_J0 = Psi * Dx_ev_J * Psi' + Tell * De * Tell';
        %----------------ce liang geng xin---xin xi xu lie-----------------
        P = DataShipStruct(i-1).Data(l,4);
        x = DataShipStruct(i-1).Data(l,1);
        y = DataShipStruct(i-1).Data(l,2);
        z = DataShipStruct(i-1).Data(l,3);
        L =  sqrt((x - X_EV_J0(1,1))*(x - X_EV_J0(1,1)) + (y - X_EV_J0(4,1))*(y - X_EV_J0(4,1)) + (z - X_EV_J0(7,1))*(z - X_EV_J0(7,1)));

        h = L + X_EV_J0(10,1);
        
        H_k = [
            -(x - X_EV_J0(1,1))/L
            0
            0
            -(y - X_EV_J0(4,1))/L
            0
            0
            -(z - X_EV_J0(7,1))/L
            0
            0
            1
            0
            0
            ]';
        
        Vz_K = P - h;
        K_k = Dx_ev_J0 * H_k' / (H_k * Dx_ev_J0 * H_k' + 4);
        X_EV_J =X_EV_J0 + K_k * Vz_K;
        Dx_ev_J = (eye(12) - K_k*H_k) * Dx_ev_J0;
        vList(i,1) = Vz_K' * 4 * Vz_K;
    end
    
    Recursive_Updating_Result(:,i) = X_EV_J;
    %DX(:,i*12-11:i*12) = Dx_ev_J;
end

%[B,L,H] = XYZ2BLH(Recursive_Updating_Result(1,:),Recursive_Updating_Result(3,:),Recursive_Updating_Result(5,:));