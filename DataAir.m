clear;clc
DataAirStruct = ReadFile('DataAir.txt');
AirCoordinate = zeros(length(DataAirStruct),4);

parfor i = 1:length(DataAirStruct)
    solut = [0 0 0 0];
    detla =[10 10 10 10];

    D = eye(DataAirStruct(i).Num) * 4;
    L = DataAirStruct(i).Data(:,4);
  
    while max(abs(detla)) > 1.0
        X = DataAirStruct(i).Data(:,1) - solut(1);
        Y = DataAirStruct(i).Data(:,2) - solut(2);
        Z = DataAirStruct(i).Data(:,3) - solut(3);
        S = sqrt(X.*X + Y.*Y + Z.*Z);
        H = [-X./S -Y./S -Z./S ones(length(X),1)];
        z = L - S - solut(4);
        detla = (H' / D * H) \ (H' / D * z);
        solut = solut + detla';
    end

    vList(i) = z' / D * z;
    AirCoordinate(i,:) = solut;
    
end
[B,L,H] = XYZ2BLH(AirCoordinate(:,1),AirCoordinate(:,2),AirCoordinate(:,3));