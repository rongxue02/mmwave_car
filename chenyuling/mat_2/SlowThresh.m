function[threshOut, reflag] = SlowThresh(dataIn, sampleN, Coeff, threshN)

%dataIn  输入数组
%sampleN  输入数组长度
%threshN  计算门限所用数据长度


%数组从小到大排列
if sampleN<0 || threshN <0
    reflag = 0;
else 
    reflag = 1;
end

DataTemp = sort(dataIn(1:sampleN));

ftemp = sum(DataTemp(1:threshN))/threshN + 100;

threshOut = ftemp*Coeff*ones(1,sampleN);