T = DataAirStruct(2).GPSTime - DataAirStruct(1).GPSTime;

Vx = (AirCoordinate(2,1) - AirCoordinate(1,1))/T;
Vy = (AirCoordinate(2,2) - AirCoordinate(1,2))/T;
Vz = (AirCoordinate(2,3) - AirCoordinate(1,3))/T;
Vt = (AirCoordinate(2,4) - AirCoordinate(1,4))/T;

Ax = ((AirCoordinate(3,1) - AirCoordinate(2,1))/T - (AirCoordinate(2,1) - AirCoordinate(1,1))/T) / T;
Ay = ((AirCoordinate(3,2) - AirCoordinate(2,2))/T - (AirCoordinate(2,2) - AirCoordinate(1,2))/T) / T;
Az = ((AirCoordinate(3,3) - AirCoordinate(2,3))/T - (AirCoordinate(2,3) - AirCoordinate(1,3))/T) / T;
At = ((AirCoordinate(3,4) - AirCoordinate(2,4))/T - (AirCoordinate(2,4) - AirCoordinate(1,4))/T) / T;

X_EV_J(1,1) = AirCoordinate(1,1) + Vx * T + 0.5*Ax*T*T;
X_EV_J(4,1) = AirCoordinate(1,2) + Vy * T + 0.5*Ay*T*T;
X_EV_J(7,1) = AirCoordinate(1,3) + Vz * T + 0.5*Az*T*T;
X_EV_J(10,1) = AirCoordinate(1,4) + Vt * T + 0.5*At*T*T;

X_EV_J(2,1) = Vx + Ax*T;
X_EV_J(5,1) = Vy + Ay*T;
X_EV_J(8,1) = Vz + Az*T;
X_EV_J(11,1) = Vt + At*T;

X_EV_J(3,1)  = Ax;
X_EV_J(6,1)  = Ay;
X_EV_J(9,1)  = Az;
X_EV_J(12,1) = At;

Dx_ev_J = eye(12) * 100;

De=diag([0.1,0.1,0.1,0.1]);

Recursive_Updating_Result = zeros(12,length(AirCoordinate));
Recursive_Updating_Result(:,1) = X_EV_J;
%DX(:,1:12) = eye(12) * 4;
vList = zeros(length(AirCoordinate),1);
T = 1;
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

for i = 2:(length(AirCoordinate)+1)
    COUNT = DataAirStruct(i-1).Num;
%     if i <= length(AirCoordinate)
%         T = DataAirStruct(i).GPSTime - DataAirStruct(i-1).GPSTime;
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
        P = DataAirStruct(i-1).Data(l,4);
        x = DataAirStruct(i-1).Data(l,1);
        y = DataAirStruct(i-1).Data(l,2);
        z = DataAirStruct(i-1).Data(l,3);
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

%[B,L,H] = XYZ2BLH(Recursive_Updating_Result(1,:),Recursive_Updating_Result(4,:),Recursive_Updating_Result(10,:));