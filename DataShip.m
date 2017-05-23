
clear;clc
DataShipStruct = ReadFile('DataShip.txt');


ShipCoordinate = zeros(length(DataShipStruct),4);
%vList = zeros(length(DataStruct),1);
parfor i = 1:length(DataShipStruct)
    solut = [0 0 0 0];
    detla =[10 10 10 10];
    %z = zeros(DataStruct(i).Num,1);
    D = eye(DataShipStruct(i).Num) * 4;
    L = DataShipStruct(i).Data(:,4);
    %H = zeros(DataStruct(i).Num,4);
    
    while max(abs(detla)) > 1.0
        X = DataShipStruct(i).Data(:,1) - solut(1);
        Y = DataShipStruct(i).Data(:,2) - solut(2);
        Z = DataShipStruct(i).Data(:,3) - solut(3);
        S = sqrt(X.*X + Y.*Y + Z.*Z);
        H = [-X./S -Y./S -Z./S ones(length(X),1)];
        z = L - S - solut(4);
        detla = (H' / D * H) \ (H' / D * z);
        solut = solut + detla';
    end
    %     X = DataStruct(i).Data(:,1) - solut(1);
    %     Y = DataStruct(i).Data(:,2) - solut(2);
    %     Z = DataStruct(i).Data(:,3) - solut(3);
    %     S = sqrt(X.*X + Y.*Y + Z.*Z);
    %     H = [-X./S -Y./S -Z./S ones(length(X),1)];
    %     z = L - S - solut(4);
    vList(i) = z' / D * z;
    ShipCoordinate(i,:) = solut;
    
end
%[B,L,H] = XYZ2BLH(ShipCoordinate(:,1),ShipCoordinate(:,2),ShipCoordinate(:,3));