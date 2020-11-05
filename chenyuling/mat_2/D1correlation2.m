function  [TotalTgtNum,Tinfo2]=D1correlation2(VlgRes,pcresult,cfar_win,SigMags,SpeedNum,RangeNum,MaxRangeCorLen)

%将VlgRes和SigMags矩阵化放在函数外处理?
%输入要求:VlgRes,SigMags均为 SpeedNum*RangeNum
%MaxRangeCorLen为最大距离相关区的长度
%要在输入之前将矩阵做好转置等操作

%对于V*R矩阵,做同一个速度门下各个距离门的相关

%struct Tinfo2
%{
%	VTargetFlag;RTargetNum;tinfo_002[]
%}
%VTargetFlag:当前速度门上目标有无标志,0:无;1:有
%RTargetNum:同一个速度门下各个距离门上目标个数
%tinfo_002[]:当前距离门上的各个目标的信息

%struct tinfo_002
%(
%	RWeight;ACCSigmags;MaxSigmags;VMax;RMax;RStart;REnd
%}
%同一个距离门上的各个目标中的一个
%RWeight:本目标的距离重心
%ACCSigmags:本目标的幅度累加
%MaxSigmags:本目标对应的脉压结果中的最大值
%VMax:本速度门号
%RMax:本目标对应的最大距离门
%RStart:距离起始单元
%REnd:距离结束单元

%Tinfo2[]
MAXTGTPERNUM=20;
%初始化Tinfo2[]
for kk1=1:SpeedNum
	Tinfo2(kk1).VTargetFlag=0;
	Tinfo2(kk1).RTargetNum=0;
	for kk2=1:MAXTGTPERNUM
		Tinfo2(kk1).tinfo_002(kk2).RWeight=0;
		Tinfo2(kk1).tinfo_002(kk2).ACCSigmags=0;
		Tinfo2(kk1).tinfo_002(kk2).MaxSigmags=0;
		Tinfo2(kk1).tinfo_002(kk2).VMax=0;
		Tinfo2(kk1).tinfo_002(kk2).RMax=0;
		Tinfo2(kk1).tinfo_002(kk2).RStart=0;
		Tinfo2(kk1).tinfo_002(kk2).REnd=0;
        Tinfo2(kk1).tinfo_002(kk2).snr=0;
	end
end% for kk1=1:RangeNum


debugCnt=0;

TotalTgtNum=0;

%考察同一个速度门下的所有距离门;
VlgResThis=zeros(1,RangeNum);
SigMagsThis=zeros(1,RangeNum);

for kkV=1:SpeedNum
	TgtNum=0;
	SigMagsThis=SigMags(kkV,:);
	VlgResThis=VlgRes(kkV,:);
	
	% [TgtNum,TgtStart,TgtEnd]=sigdetect_m_n(VlgRes,Threshold,SigLen)
	%当前速度门下所有距离门的检测
	[TgtNum,TgtStart,TgtEnd] = sigdetect_m_n(VlgResThis,0.5,RangeNum);
	%目标个数限制
    if(TgtNum > MAXTGTPERNUM)
        TgtNum = MAXTGTPERNUM;
    end
    
	if TgtNum>0
		%fvar1:幅度加权累加值
		%fvar2:幅度直接累加值
		fvar1=0;
		fvar2=0;
        mmTgt=0;
		SigMaxValueArray=zeros(1,TgtNum);
		%目标个数
		for kkTgt=1:TgtNum
			fvar1=0;fvar2=0;
			%每个目标所占的距离门		
			for kkTgtIndex=TgtStart(kkTgt):TgtEnd(kkTgt)
				fvar1=fvar1+kkTgtIndex*SigMagsThis(kkTgtIndex);
				fvar2=fvar2+SigMagsThis(kkTgtIndex);
            end

            %debugCnt=debugCnt+1
            [SigMaxValue,SigMaxIndex]=max(SigMagsThis(TgtStart(kkTgt):TgtEnd(kkTgt)));
			SigMaxValueArray(kkTgt)=SigMaxValue;
			
            if 1%(TgtEnd(kkTgt)-TgtStart(kkTgt))<=MaxRangeCorLen
                mmTgt = mmTgt+1;
                Tinfo2(kkV).tinfo_002(mmTgt).RWeight=fvar1/fvar2;
                Tinfo2(kkV).tinfo_002(mmTgt).ACCSigmags=fvar2;
                Tinfo2(kkV).tinfo_002(mmTgt).MaxSigmags=SigMaxValue;
                Tinfo2(kkV).tinfo_002(mmTgt).VMax=kkV;
                Tinfo2(kkV).tinfo_002(mmTgt).RMax=TgtStart(kkTgt)+SigMaxIndex-1;
                Tinfo2(kkV).tinfo_002(mmTgt).snr=SigMags(kkV,(TgtStart(kkTgt)+SigMaxIndex-1))/pcresult(kkV,(TgtStart(kkTgt)+SigMaxIndex-1))*cfar_win;
                Tinfo2(kkV).tinfo_002(mmTgt).RStart=TgtStart(kkTgt);
                Tinfo2(kkV).tinfo_002(mmTgt).REnd=TgtEnd(kkTgt);
            end;
% 			Tinfo2(kkV).tinfo_002(kkTgt).RWeight=fvar1/fvar2;
% 			Tinfo2(kkV).tinfo_002(kkTgt).ACCSigmags=fvar2;
% 			Tinfo2(kkV).tinfo_002(kkTgt).MaxSigmags=SigMaxValue;
% 			Tinfo2(kkV).tinfo_002(kkTgt).VMax=kkV;
% 			%Tinfo2(kkV).tinfo_002(kkTgt).RMax=SigMaxIndex;
%             Tinfo2(kkV).tinfo_002(kkTgt).RMax=TgtStart(kkTgt)+SigMaxIndex-1;
% 			Tinfo2(kkV).tinfo_002(kkTgt).RStart=TgtStart(kkTgt);
% 			Tinfo2(kkV).tinfo_002(kkTgt).REnd=TgtEnd(kkTgt);
        end
		
        if mmTgt>0
            Tinfo2(kkV).VTargetFlag=1;
        else
            Tinfo2(kkV).VTargetFlag=0;
        end
		Tinfo2(kkV).RTargetNum=mmTgt;
		TotalTgtNum=TotalTgtNum+mmTgt;
		
        %%%%%%%%%%%%%%%%%距离维不需要进行首尾合成处理 林家豪 20190731%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		%将由于FFT引起的首尾的分裂合成(liu);首尾边界处理
% 		%同一个速度门上不止一个目标
% 		%对第一个和最后一个目标进行特殊处理
% 		if mmTgt>1
% 			TempTgtStart=Tinfo2(kkV).tinfo_002(1).RStart;
% 			TempTgtEnd=Tinfo2(kkV).tinfo_002(mmTgt).REnd;
% 			
% 			%存在首尾分裂的情况
% 			if TempTgtStart==1 && TempTgtEnd==RangeNum
% 				%bandfvar1:幅度加权累加值
% 				%bandfvar2:幅度直接累加值
% 				bandfvar1=0;
% 				bandfvar2=0; 
% 				SigMaxValue=0;
% 				SigMaxIndex=0;
% 				TempNewTgtStart=Tinfo2(kkV).tinfo_002(mmTgt).RStart-RangeNum;
% 				TempNewTgtEnd=Tinfo2(kkV).tinfo_002(1).REnd;
% 				
% 				for kkTgtIndex=TempNewTgtStart:TempNewTgtEnd
% 					NewIndex=mod(kkTgtIndex+RangeNum,RangeNum);
% 					bandfvar1=bandfvar1+kkTgtIndex*SigMagsThis(NewIndex);
% 					bandfvar2=SigMagsThis(NewIndex);
% 					if SigMagsThis(NewIndex)>SigMaxValue
% 						SigMaxValue=SigMagsThis(NewIndex);
% 						SigMaxIndex=kkTgtIndex;
% 					end 	
% 				end
% 				%for
% 				
% 				SigMaxValueArray(1)=SigMaxValue;
% 				%该目标的重心
% 				TempRWeight=bandfvar1/bandfvar2;
% 				
% 				if TempRWeight<0
% 					Tinfo2(kkV).tinfo_002(1).RWeight=TempRWeight+RangeNum;
% 				else
% 					Tinfo2(kkV).tinfo_002(1).RWeight=TempRWeight;
% 				end
% 				
% 				Tinfo2(kkV).tinfo_002(1).ACCSigmags=bandfvar2;
% 				Tinfo2(kkV).tinfo_002(1).MaxSigmags=SigMaxValue;
% 				Tinfo2(kkV).tinfo_002(1).VMax=kkV;
% 				Tinfo2(kkV).tinfo_002(1).RMax=mod(SigMaxIndex+RangeNum,RangeNum);
% 				Tinfo2(kkV).tinfo_002(1).RStart=mod(TempNewTgtStart+RangeNum,RangeNum);
% 				Tinfo2(kkV).tinfo_002(1).REnd=TempNewTgtEnd;
% 				Tinfo2(kkV).VTargetFlag=1;
% 				%合并掉一个目标
% 				Tinfo2(kkV).RTargetNum=mmTgt-1;
% 				mmTgt=mmTgt-1;
% 				TotalTgtNum=TotalTgtNum-1;
% 
% 			end
% 			%TempTgtStart==1 && TempTgtEnd==RangeNum	
% 		end
		%TgtNum>1
	else
		%没发现目标
		Tinfo2(kkV).VTargetFlag=0;
		
	end
	%TgtNum>0
	if 0 
	%----------------------------
	%剔除滤波器旁瓣输出的虚假目标
	% if Tinfo2(kkV).RTargetNum>0
	if TgtNum>0
		SigMaxValueThis=max(SigMaxValueArray);
		%将旁瓣高度设为最大值的0.1%20*LOG0.1=-20dB
		SigSideMag=0.1*SigMaxValueThis;
		
		NewIndex = 1;
		NewIndexRef=0;
		
		for kkTgt=1:Tinfo2(kkV).RTargetNum
			%去掉小旁瓣
			if Tinfo2(kkV).tinfo_002(kkTgt).MaxSigmags>SigSideMag
				Tinfo2(kkV).tinfo_002(NewIndex)=Tinfo2(kkV).tinfo_002(kkTgt);
				NewIndex=NewIndex+1;
				NewIndexRef=NewIndexRef+1;
			end
		end
		%更新总的目标数
		TotalTgtNum=TotalTgtNum-Tinfo2(kkV).RTargetNum+NewIndexRef;
		%更新本距离门的目标个数
		Tinfo2(kkV).RTargetNum=NewIndexRef;
	end
	%Tinfo2(kkV).RTargetNum>0
    end
end
%kkV=1:SpeedNum

end