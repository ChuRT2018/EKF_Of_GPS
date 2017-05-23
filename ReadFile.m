function [ Data ] = ReadFile( file )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [name,r1,r2,r3,r4] = textread(file,'%s %f %f %f %f','commentstyle','matlab');
    RawData = [r1 r2 r3 r4];
     
    Length = 0;
    
    Last.GPSWeek = 0;
    Last.GPSTime = 0;
    % use struct 'DataStatic'
    for i = 1 : length(name)
        if isequal(cell2mat(name(i,1)),'#')
            Length = Length + 1;
            
            if r2(i) == Last.GPSWeek && r3(i) == Last.GPSTime
                continue;
            end
            
            m = 1;
            while ~isequal(cell2mat(name(i+m,1)),'#')
                if (i + m) == length(name)
                    break;
                end
                m = m + 1;
            end
            
            Data(Length).ID = r1(i);
            Data(Length).GPSWeek = r2(i);
            Data(Length).GPSTime = r3(i);
            Data(Length).Num = m - 1;
            Data(Length).Data = zeros(m - 1,4);
            % Storing data
            for l = 1: Data(Length).Num
                Data(Length).Data(l,:) = RawData(i+l,:);
            end
            Last = Data(Length);
        end   
    end
end

