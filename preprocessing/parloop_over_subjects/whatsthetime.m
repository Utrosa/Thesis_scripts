function time = whatsthetime
T = datetime('now');
if length(num2str(hour(T))) == 1
    hr = ['0',num2str(hour(T))];
else 
    hr = num2str(hour(T));
end
if length(num2str(minute(T))) == 1
    min = ['0',num2str(minute(T))];
else
    min = num2str(minute(T));
end
time = [hr,':',min];
end