%%  ReadWindData  by yy 2021.8.27
function [data_Wind] = ReadWindData(data,Filedir)

wind_num = data.Wind.wind_Number; %number of wind units

wind = xlsread(Filedir); %read wind data
data = wind(1:end,2:end);
data = data(:,9:8+wind_num);
data = rmmissing(data,1);
data = data * 0.0009;
data(data < 1e-3) = 0;

for i = 1:wind_num
    data_Wind{i} = reshape(data(:,i),24,[]);
end

end

