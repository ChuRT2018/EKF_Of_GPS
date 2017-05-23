[name,r1,r2,r3,r4] = textread('DataStatic.txt','%s %f %f %f %f','commentstyle','matlab');
RawData = [r1 r2 r3 r4];

global DataStatic;
Length = 0;

% use static 'DataStatic'
for i = 1 : length(name)
    if isequal(cell2mat(name(i,1)),'#')
        Length = Length + 1;
        m = 1;
        while ~isequal(cell2mat(name(i+m,1)),'#')
            if i + m == length(name)
                break;
            end
            m = m + 1;
        end
        DataStatic(Length).ID = r1(i);
        DataStatic(Length).GPSWeek = r2(i);
        DataStatic(Length).GPSTime = r3(i);
        DataStatic(Length).Num = m - 1;
        DataStatic(Length).Data = zeros(r4(i),4);
        % Storing data
        for l = 1: DataStatic(Length).Num
            DataStatic(Length).Data(l,:) = RawData(i+l,:);
        end
    end
end



