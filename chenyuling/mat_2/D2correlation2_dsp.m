function [Plot,CurrentTarNum]=D2correlation2_dsp(Tinfo2,SpeedNum,MaxPlotNum,varVCONDITION)
%做二维相关;Tinfo为对应的一维相关的结果;SpeedNum*RangeNum

%struct Plot
%{
%	V;VStart;VEnd;
%	Sigmags;RMax;VMax;
%	R;RStart;REnd;
%	Azimuth;Time;Flag
%}
%V:速度重心;VStart:速度起始单元;VEnd:速度结束单元;
%Sigmags:信号幅度最大值;RMax:信号幅度最大值所在的距离门;VMax:信号幅度最大值所在的速度门;
%R:距离重心;RStart:距离起始单元;REnd:距离结束单元;
%Azimuth:当前帧的方位角;Time:当前帧的时间;Flag:相关标志
%varVCONDITION:为最大速度相关区的长度

%需要尤其注意的是:我们的输入结果已经是一维相关的结果了,已经脱离了点的概念,在我们面前的是一个个的目标


%MaxPlotNum:最多可以检测出多少个目标(input)
varMAXPLOTNUM=MaxPlotNum;

%初始化结果区
for kk1=1:MaxPlotNum
	Plot(kk1).V=0;
	Plot(kk1).VStart=0;
	Plot(kk1).VEnd=0;
	
    Plot(kk1).snr=0;
	Plot(kk1).Sigmags=0;
	Plot(kk1).RMax=0;
	Plot(kk1).VMax=0;
	
	Plot(kk1).R=0;
	Plot(kk1).RStart=0;
	Plot(kk1).REnd=0;
	
    Plot(kk1).ACCSigmags = 0;
% 	Plot(kk1).Azimuth=0;
	Plot(kk1).Time=0;
	Plot(kk1).Flag=0;
end
CurrentTarNum=0;
%将Tinfo2结构体中的目标个数单独取出来，单独赋值到数组中
cellTemp = {Tinfo2.RTargetNum};
Index_1D = cell2mat(cellTemp);
%找出所有目标的
Index_R = find(Index_1D > 0);
%能检测到目标的速度滤波器个数
ftemp1 = length(Index_R);


%--------Plot中按照第一个有目标的速度滤波器填充数据
if ftemp1 == 0
    %一维检测没有检测到目标
    CurrentTarNum =0;
else
    Valid1st = Index_R(1); 
    % Valid1st为当前的速度滤波器号
    for kk2 = 1:Index_1D(Valid1st)
        if CurrentTarNum >= varMAXPLOTNUM
            break;
        end	
        CurrentTarNum = CurrentTarNum + 1;
        Plot(CurrentTarNum).V         = Tinfo2(Valid1st).tinfo_002(kk2).VMax;
        Plot(CurrentTarNum).VStart    = Tinfo2(Valid1st).tinfo_002(kk2).VMax;
        Plot(CurrentTarNum).VEnd      = Tinfo2(Valid1st).tinfo_002(kk2).VMax;

        Plot(CurrentTarNum).snr       = Tinfo2(Valid1st).tinfo_002(kk2).snr;
        Plot(CurrentTarNum).Sigmags   = Tinfo2(Valid1st).tinfo_002(kk2).MaxSigmags;
        Plot(CurrentTarNum).RMax      = Tinfo2(Valid1st).tinfo_002(kk2).RMax;
        Plot(CurrentTarNum).VMax      = Tinfo2(Valid1st).tinfo_002(kk2).VMax;

        Plot(CurrentTarNum).R         = Tinfo2(Valid1st).tinfo_002(kk2).RWeight;
        Plot(CurrentTarNum).RStart    = Tinfo2(Valid1st).tinfo_002(kk2).RStart;
        Plot(CurrentTarNum).REnd      = Tinfo2(Valid1st).tinfo_002(kk2).REnd;

        Plot(CurrentTarNum).ACCSigmags= Tinfo2(Valid1st).tinfo_002(kk2).ACCSigmags;
    end


    for ii1=2:ftemp1
        kk1 = Index_R(ii1); 
        % kk1为当前的速度滤波器号
        % 如果二维点迹数为零，则用一维点迹初始化

        %对每一个一维点迹，用已存在的二维点迹进行相关，相关不上的建立新的二维点迹
        flagA = ones(1,varMAXPLOTNUM);%点迹相关上标志
        for kk2=1:Tinfo2(kk1).RTargetNum    % kk2为该速度滤波器上一维点迹序号
            flg = 0;
            a2 = Tinfo2(kk1).tinfo_002(kk2).RStart;
            b2 = Tinfo2(kk1).tinfo_002(kk2).REnd;
            for kk3 = 1:CurrentTarNum       % CurrentTarNum为当前二维点迹
                a1 = Plot(kk3).RStart;
                b1 = Plot(kk3).REnd;
                if flagA(kk3) == 0
                    continue;
                end
                % 判断速度滤波器是否相邻
                if (kk1-Plot(kk3).VEnd) == 1
                    % 距离交叠判断
                    if ((a2>=a1)&&(a2<=b1)) || ((a1>=a2)&&(a1<=b2))
                        % 相关成功
                        flg = 1; 
                        Plot(kk3).RStart    = min(a1, a2);
                        Plot(kk3).REnd      = max(b1, b2);
                        if Plot(kk3).Sigmags < Tinfo2(kk1).tinfo_002(kk2).MaxSigmags
                            Plot(kk3).Sigmags = Tinfo2(kk1).tinfo_002(kk2).MaxSigmags;
                            Plot(kk3).RMax    = Tinfo2(kk1).tinfo_002(kk2).RMax;
                            Plot(kk3).VMax    = kk1; 
                        end
                        RAcc                 = Plot(kk3).R*Plot(kk3).ACCSigmags + Tinfo2(kk1).tinfo_002(kk2).RWeight*Tinfo2(kk1).tinfo_002(kk2).ACCSigmags;
                        VAcc                  = Plot(kk3).V*Plot(kk3).ACCSigmags + Tinfo2(kk1).tinfo_002(kk2).ACCSigmags*kk1;
                        Plot(kk3).ACCSigmags  = Plot(kk3).ACCSigmags + Tinfo2(kk1).tinfo_002(kk2).ACCSigmags;
                        Plot(kk3).R           = RAcc / Plot(kk3).ACCSigmags;
                        Plot(kk3).V           = VAcc / Plot(kk3).ACCSigmags;
                        Plot(kk3).VEnd        = kk1;
                        Plot(kk3).snr = max(Plot(kk3).snr,Tinfo2(kk1).tinfo_002(kk2).snr);
                        break;
                    end
                else
                    flagA(kk3) = 0;          % 若不相邻，则将该二维点迹相关标志位清零，这样对于当前滤波器上的所有一维点迹都不必再与此二维点迹做相关
                end
            end
            % 若没有相关上，则建立新点迹
            if flg == 0
                if CurrentTarNum >= varMAXPLOTNUM
                    break;
                end	
                CurrentTarNum = CurrentTarNum + 1;
                Plot(CurrentTarNum).V         = Tinfo2(kk1).tinfo_002(kk2).VMax;
                Plot(CurrentTarNum).VStart    = Tinfo2(kk1).tinfo_002(kk2).VMax;
                Plot(CurrentTarNum).VEnd      = Tinfo2(kk1).tinfo_002(kk2).VMax;

                Plot(CurrentTarNum).snr   = Tinfo2(kk1).tinfo_002(kk2).snr;
                Plot(CurrentTarNum).Sigmags   = Tinfo2(kk1).tinfo_002(kk2).MaxSigmags;
                Plot(CurrentTarNum).RMax      = Tinfo2(kk1).tinfo_002(kk2).RMax;
                Plot(CurrentTarNum).VMax      = Tinfo2(kk1).tinfo_002(kk2).VMax;

                Plot(CurrentTarNum).R         = Tinfo2(kk1).tinfo_002(kk2).RWeight;
                Plot(CurrentTarNum).RStart    = Tinfo2(kk1).tinfo_002(kk2).RStart;
                Plot(CurrentTarNum).REnd      = Tinfo2(kk1).tinfo_002(kk2).REnd;

                Plot(CurrentTarNum).ACCSigmags= Tinfo2(kk1).tinfo_002(kk2).ACCSigmags;
            end
        end
    end
    % for ii = 1:CurrentTarNum
    %         aq(ii) = Plot(ii).Sigmags;
    %     end
    %     [ak,nk] = sort(aq,'descend');
    %     Plot1 = Plot(nk);
    %%%%%%%%%%%%%%速度边界处理%%%%%%%%%%%%%%%
    %前面提出来的点迹按多普勒维从小到大排列
    flagB = ones(1,varMAXPLOTNUM);%点迹有效标志
    for kki=1:CurrentTarNum-1
        if (flagB(kki) == 1)&&(Plot(kki).VStart == 1)%该点迹有效
            for kkj=(kki+1):CurrentTarNum
                if (flagB(kkj) == 1)&&(Plot(kkj).VEnd == SpeedNum)&&...
                   (((Plot(kki).RStart>= Plot(kkj).RStart)&&(Plot(kki).RStart<=Plot(kkj).REnd))||(( Plot(kkj).RStart>=Plot(kki).RStart)&&( Plot(kkj).RStart<=Plot(kki).REnd)))  
                        %对于这两个目标求一个新的距离重心
                        fweigth0=Plot(kki).R*Plot(kki).ACCSigmags;
                        fweigth1=Plot(kkj).R*Plot(kkj).ACCSigmags;
                        Plot(kki).R=(fweigth0+fweigth1)/(Plot(kki).ACCSigmags+Plot(kkj).ACCSigmags);
                        %求二维相关区的距离维的最小值与最大值
                        if Plot(kkj).RStart<Plot(kki).RStart
                            Plot(kki).RStart=Plot(kkj).RStart;
                        end
                        if Plot(kkj).REnd>Plot(kki).REnd
                            Plot(kki).REnd=Plot(kkj).REnd;
                        end
                        %得到幅度最大值及坐标、取最大范围得到速度起始和结束坐标
                        if Plot(kkj).Sigmags>=Plot(kki).Sigmags
                            Plot(kki).Sigmags = Plot(kkj).Sigmags;
                            Plot(kki).RMax = Plot(kkj).RMax;
                            Plot(kki).VMax = Plot(kkj).VMax;                               
                        end
                        %kki在左,kkj在右
                        %对于这两个目标求一个新的速度重心
                        fweigth0=Plot(kki).V*Plot(kki).ACCSigmags;
                        fweigth1=( Plot(kkj).V-SpeedNum )*Plot(kkj).ACCSigmags;
                        Plot(kki).V=(fweigth0+fweigth1)/(Plot(kki).ACCSigmags+Plot(kkj).ACCSigmags);
                        Plot(kki).V = mod((Plot(kki).V+SpeedNum),SpeedNum);
                        %求二维相关区的速度维的最小值与最大值
                        Plot(kki).VStart = Plot(kkj).VStart;
                        Plot(kki).VEnd = Plot(kki).VEnd;
                        %对于这两个目标求一个新的幅度累加和,在计算完重心后,否则会错误
                        Plot(kki).ACCSigmags = Plot(kki).ACCSigmags + Plot(kkj).ACCSigmags;
                        Plot(kki).snr = max(Plot(kki).snr,Plot(kkj).snr);
                        flagB(kkj)=0;%被融合点迹设置标志无效
                end
            end%for kkj=(kki+1):SpeedNum
        end%(flagB(kki) == 1)&&(Plot(kki).VStart == 1)%该点迹有效
    end%for kki=1:SpeedNum-1

    kksave = 0;
    for kki=1:CurrentTarNum
        if (Plot(kki).VStart <= Plot(kki).VEnd)
            vlen = abs(Plot(kki).VEnd-Plot(kki).VStart+1);
        else
            vlen = abs(Plot(kki).VEnd-Plot(kki).VStart+SpeedNum+1);
        end
        if (flagB(kki) == 1)&&(vlen<=varVCONDITION)
            kksave = kksave+1;%有效点迹数
            Plot(kksave) = Plot(kki);
        end
    end
    for kki=kksave+1:CurrentTarNum
        Plot(kki).V=0;
        Plot(kki).VStart=0;
        Plot(kki).VEnd=0;

        Plot(kki).snr=0;
        Plot(kki).Sigmags=0;
        Plot(kki).RMax=0;
        Plot(kki).VMax=0;

        Plot(kki).R=0;
        Plot(kki).RStart=0;
        Plot(kki).REnd=0;

        Plot(kki).ACCSigmags = 0;
    % 	Plot(kki).Azimuth=0;
        Plot(kki).Time=0;
    end
    CurrentTarNum = kksave;

end  %if ftemp1 > 0
end