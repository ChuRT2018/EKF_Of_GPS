T = DataShipStruct(2).GPSTime - DataShipStruct(1).GPSTime;

Vx = (ShipCoordinate(2,1) - ShipCoordinate(1,1))/T;
Vy = (ShipCoordinate(2,2) - ShipCoordinate(1,2))/T;
Vz = (ShipCoordinate(2,3) - ShipCoordinate(1,3))/T;
Vt = (ShipCoordinate(2,4) - ShipCoordinate(1,4))/T;

X_EV_J(1,1) = ShipCoordinate(1,1) + Vx * T;
X_EV_J(3,1) = ShipCoordinate(1,2) + Vy * T;
X_EV_J(5,1) = ShipCoordinate(1,3) + Vz * T;
X_EV_J(7,1) = ShipCoordinate(1,4) + Vt * T;

X_EV_J(2,1) = Vx;
X_EV_J(4,1) = Vy;
X_EV_J(6,1) = Vz;
X_EV_J(8,1) = Vt;

Dx_ev_J = eye(8) * 100;

De=diag([0.01,0.01,0.01,0.01]);

Recursive_Updating_Result = zeros(8,length(ShipCoordinate));
Recursive_Updating_Result(:,1) = X_EV_J;
DX(:,1:8) = eye(8) * 4;
vList = zeros(length(ShipCoordinate),1);

for i = 2:(length(ShipCoordinate)+1)
    COUNT = DataShipStruct(i-1).Num;
    if i <= length(ShipCoordinate)
        T = DataShipStruct(i).GPSTime - DataShipStruct(i-1).GPSTime;
    end
    
    Psi = [
        1       T      0      0        0     0       0    0
        0       1      0      0        0     0       0    0
        0       0      1      T        0     0       0    0
        0       0      0      1        0     0       0    0
        0       0      0      0        1     T       0    0
        0       0      0      0        0     1       0    0
        0       0      0      0        0     0       1    T
        0       0      0      0        0     0       0    1
        ];
    
    Tell = [
        T 0 0 0 
        1 0 0 0
        0 T 0 0
        0 1 0 0
        0 0 T 0
        0 0 1 0
        0 0 0 T
        0 0 0 1
        ];
    
    for l = 1:COUNT
        X_EV_J0 =  Psi * X_EV_J;
        Dx_ev_J0 = Psi * Dx_ev_J * Psi' + Tell * De * Tell';
        %----------------ce liang geng xin---xin xi xu lie-----------------
        P = DataShipStruct(i-1).Data(l,4);
        x = DataShipStruct(i-1).Data(l,1);
        y = DataShipStruct(i-1).Data(l,2);
        z = DataShipStruct(i-1).Data(l,3);
        
        h = sqrt((x - X_EV_J0(1,1))*(x - X_EV_J0(1,1)) + (y - X_EV_J0(3,1))*(y - X_EV_J0(3,1)) + (z - X_EV_J0(5,1))*(z - X_EV_J0(5,1))) + X_EV_J0(7,1);
        
        H_k = [
            -(x - X_EV_J0(1,1))/sqrt((x - X_EV_J0(1,1))*(x - X_EV_J0(1,1)) + (y - X_EV_J0(3,1))*(y - X_EV_J0(3,1)) + (z - X_EV_J0(5,1))*(z - X_EV_J0(5,1)))
            0
            -(y - X_EV_J0(3,1))/sqrt((x - X_EV_J0(1,1))*(x - X_EV_J0(1,1)) + (y - X_EV_J0(3,1))*(y - X_EV_J0(3,1)) + (z - X_EV_J0(5,1))*(z - X_EV_J0(5,1)))
            0
            -(z - X_EV_J0(5,1))/sqrt((x - X_EV_J0(1,1))*(x - X_EV_J0(1,1)) + (y - X_EV_J0(3,1))*(y - X_EV_J0(3,1)) + (z - X_EV_J0(5,1))*(z - X_EV_J0(5,1)))
            0
            1
            0
            ]';
        
        Vz_K = P - h;
        K_k = Dx_ev_J0 * H_k' / (H_k * Dx_ev_J0 * H_k' + 4);
        X_EV_J =X_EV_J0 + K_k * Vz_K;
        Dx_ev_J = (eye(8) - K_k*H_k) * Dx_ev_J0;
        vList(i,1) = Vz_K' * 4 * Vz_K;
    end
    
    Recursive_Updating_Result(:,i) = X_EV_J;
    DX(:,i*8-7:i*8) = Dx_ev_J;
end

%[B,L,H] = XYZ2BLH(Recursive_Updating_Result(1,:),Recursive_Updating_Result(3,:),Recursive_Updating_Result(5,:));