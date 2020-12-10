function [binnedvector] = BinData(Vector,Bins,Min,Max)
%first argument = vector to be binned, 
% second arg = number of bins
% third and fourth args = limits
%returns a vector containing the bin indices corresponding to each element
%in Vector
%%
vector = Vector;
%Max = nanmax(vector);
%Min = nanmin(vector);
Range = Max-Min;
BinWidth = Range/Bins;
binnedvector =[]; %vector to fill with binning info
for value = 1:length(vector)
    datavalue = vector(value);
    distance = datavalue - Min;
    if distance == 0; %if the value = the minimum value, then it's in bin 1
        valuebin = 1;
    elseif datavalue == NaN
        valuebin = NaN;
    else
        valuebin = ceil(distance/BinWidth);
    end
    binnedvector(value) = valuebin; %fill the output vector
end

