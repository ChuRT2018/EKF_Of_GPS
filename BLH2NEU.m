%1. Z rotation L degree and now OXZ is coinside with the LH
%2. Y rotation 90-B degree and Z pass though point P
%3. move the origin to the point P
%4. reverse X and get EUK
NEU = zeros(3,length(L));
X0 = -1974855.2452;
Y0 = 4589660.2291;
Z0 = 3953376.2647;
for i = 1 : length(Coordinate)
    Wz = [
        cos(L(i))      sin(L(i))  0
        -sin(L(i))     cos(L(i))  0
           0               0        1
        ];
    
    
    
    Wy = [
        cos(90-B(i))    0   -sin(90-B(i))
            0           1       0
        sin(90-B(i))    0   cos(90-B(i))
        ];
    
    
    NEU(:,i) = Wy * Wz * [Coordinate(i,1) - X0;Coordinate(i,2) - Y0;Coordinate(i,3) - Z0];
end

NEU(1,:) = -NEU(1,:);