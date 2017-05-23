clear;clc;
DataStruct = ReadFile('DataStatic.txt');


Coordinate = zeros(length(DataStruct),4);
vList = zeros(length(DataStruct),1);
parfor i = 1:length(DataStruct)
    solut = [0 0 0 0];
    detla =[10 10 10 10];
    z = zeros(DataStruct(i).Num,1);
    D = eye(DataStruct(i).Num) * 4;
    L = DataStruct(i).Data(:,4);
    H = zeros(DataStruct(i).Num,4);
    
    while max(abs(detla)) > 1.0
        X = DataStruct(i).Data(:,1) - solut(1);
        Y = DataStruct(i).Data(:,2) - solut(2);
        Z = DataStruct(i).Data(:,3) - solut(3);        
        S = sqrt(X.*X + Y.*Y + Z.*Z);
        H = [-X./S -Y./S -Z./S ones(length(X),1)];
        z = L - S - solut(4);
        detla = (H' / D * H) \ (H' / D * z);
        solut = solut + detla';       
    end
    
    X = DataStruct(i).Data(:,1) - solut(1);
    Y = DataStruct(i).Data(:,2) - solut(2);
    Z = DataStruct(i).Data(:,3) - solut(3);        
    S = sqrt(X.*X + Y.*Y + Z.*Z);
    H = [-X./S -Y./S -Z./S ones(length(X),1)];
    z = L - S - solut(4);
    vList(i) = z' / D * z;
    Coordinate(i,:) = solut;

end
[B,L,H] = XYZ2BLH(Coordinate(:,1),Coordinate(:,2),Coordinate(:,3));