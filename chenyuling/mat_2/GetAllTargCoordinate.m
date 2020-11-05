function [Rbuf,Vbuf,Ampbuf]=GetAllTargCoordinate(flag,Plot)
%取出二维相关后各目标坐标信息
%flag=1使用最大值,=0 使用重心
%struct Plot
		%{
		%	V;VStart;VEnd;
		%	Sigmags;RMax;VMax;
		%	R;RStart;REnd;
		%	Azimuth;Time;Flag
		%}
%取出所有点距离坐标
if flag == 1
    %使用最大值
    celltemp1 = {Plot(:).RMax};
    datatemp1 = cell2mat(celltemp1);
else
    %使用重心
    celltemp0 = {Plot(:).R};
    %注意：中心计算得到的是小数，在后续计算中需要用到点迹坐标，因此需要将小数转成整数，选择左侧作为坐标
    datatemp1 = floor(cell2mat(celltemp0));

end
Rbuf = datatemp1;

%取出二维相关后各目标速度坐标
celltemp2 = {Plot(:).VMax};
datatemp2 = cell2mat(celltemp2);
Vbuf = datatemp2;

%各个目标幅度
celltemp3 = {Plot(:).Sigmags};
Ampbuf = cell2mat(celltemp3);
