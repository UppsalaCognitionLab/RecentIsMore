function Y = get_included_subjects(expnr)
if expnr==1
    Y = [1:23 25:30];  % one subject excluded (see paper for details)
elseif expnr==2
    Y = [1:32 34:56 58:63 65:85]; % two subjects excluded (see paper for details)
elseif expnr==3
    Y = 1:20;
end