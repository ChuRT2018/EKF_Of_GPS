% clear;clc;
% [name,r1,r2,r3,r4] = textread('DataStatic1.txt','%s %f %f %f %f','commentstyle','matlab');
% RawData = [r1 r2 r3 r4];
% 
% global Data;
% Length = 0;
% 
% Last.GPSWeek = 0;
% Last.GPSTime = 0;
% % use struct 'DataStatic'
% for i = 1 : length(name)
%     if isequal(cell2mat(name(i,1)),'#')
%         if r1(i) <= 15 || (r1(i) <= 11900 && r1(i) >= 11800) || (r1(i) <= 27500 && r1(i) >= 27600) 
%             continue;
%         end
%         
%         Length = Length + 1;
%         
%         if r2(i) == Last.GPSWeek && r3(i) == Last.GPSTime
%             continue;
%         end
%         
%         m = 1;
%         while ~isequal(cell2mat(name(i+m,1)),'#')
%             if i + m == length(name)
%                 break;
%             end
%             m = m + 1;
%         end
%         
%         Data(Length).ID = r1(i);
%         Data(Length).GPSWeek = r2(i);
%         Data(Length).GPSTime = r3(i);
%         Data(Length).Num = m - 1;
%         Data(Length).Data = zeros(m - 1,4);
%         % Storing data
%         for l = 1: Data(Length).Num
%             Data(Length).Data(l,:) = RawData(i+l,:);
%         end
%         Last = Data(Length);
%     end
% end
% 
% DataStruct = Data;
clear;clc
DataStruct = ReadFile('DataStatic1.txt');
Coordinate = zeros(length(DataStruct),4);
vList = zeros(length(DataStruct),1);
parfor i = 1:length(DataStruct)
    solut = [0 0 0 0];
    detla =[10 10 10 10];
    %z = zeros(DataStruct(i).Num,1);
    D = eye(DataStruct(i).Num) * 4;
    L = DataStruct(i).Data(:,4);
    %H = zeros(DataStruct(i).Num,4);
    
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
%XYZ2BLH;