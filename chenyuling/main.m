clear all;close all;clc;

%检测最大目标数量
TargNumMax = 64;

%包含文件，函数，定义结构体
addpath('SubFunction');

%% paramRadar

lightspeed                      = 299792458;                       %光速
%mmwave_studio配置参数
txNum                           = int32(2);                            %发射天线
rxNum                           = int32(4);                            %接收天线
sampleRate                    = 10e6;                                %采样率
rangeSmpNum              = int32(128);                          %采样点数
chirpNum                      = int32(128);                         %chirp个数
rangeSmpTime              = double(rangeSmpNum)/sampleRate;  %采样时间
centerFreq                     = 77e9;                              %中心频率
freqSlope                      = 29.982e12;                         %调频斜率，单位Hz/s
idleTime                        = 30e-6;                            %空闲时间
ADCstartTime                = 6e-6;                            %ADC开始时间
rampEndTime                = 40e-6;                           %调频持续时间
framePeriod                  = 40e-3;                           %驻留周期
FrameNum                    = 1000;                             %采集驻留个数

%通过配置参数计算得到

Rx_n                         = rxNum * txNum;           %等效虚拟天线数目

lamda                        = lightspeed/centerFreq;   %波长
bandwidth                    = freqSlope*rampEndTime;	%带宽
chirpRepTime                 = double(txNum)*(idleTime + rampEndTime);	%一个chirp的时间
chirpRepFreq                 = 1/(chirpRepTime);	%FMCW波形频率
chirpTime                    = double(chirpNum)*chirpRepTime;	%所有chirp占用的时间
rangeRes                     = 299792458/2/freqSlope*sampleRate/double(rangeSmpNum);    %距离分辨率，c/2/miu*Fs/SampleN=c/2B
velocityRes                  = 299792458/centerFreq/2/chirpRepTime/double(chirpNum); %速度分辨率，c/fc/2/ChirpTime/ChirpN
maxRange                     = 299792458/2/freqSlope*sampleRate;		%测距最大值，c/2/miu*Fs
maxVelocity                  = 299792458/centerFreq/4*chirpRepFreq;    	%测速最大值c/fc/4/ChirpTime
fIFMax                       = freqSlope*2*maxRange/lightspeed;
WorkTime                     = framePeriod - chirpTime;     %chirp后剩余时间，用来处理

%根据1642检测角度范围120度推算出接收天线间距
distanceRX                   = lamda/2;     %接收天线间距
distanceTX                   = distanceRX*4;     %发射天线间距

SampleM = double(rangeSmpNum);
ChirpN = double(chirpNum);
RxP = double(rxNum);
dPerRx = distanceRX;
%3维FFT点数设置
RFFTNum =1* power(2,ceil(log2(SampleM)));       %距离维FFT点数
MTDFftNum =1* power(2,ceil(log2(ChirpN)));      %多普勒点数
FFTNumdbf =128;                              %不同接收通道做FFT点数，求角度
%每个点对应的信息            
rangeCell               = maxRange/RFFTNum;	%FFT后chirp内每个点代表的距离
velocityCell            = 2*maxVelocity/MTDFftNum;	%FFT后chirp内每个点代表的速度,相位范围[-1,1],因此计算步进的时候要*2
angleCell = 2*asin(lamda/(2*dPerRx)) *180/pi/FFTNumdbf;    %角度坐标%相位范围[-1,1],因此计算步进的时候要*2
%坐标
X_All = ( -RFFTNum/2 )*rangeCell:rangeCell:(RFFTNum/2-1)*rangeCell;    %X轴坐标,有负数
%由于距离没有负值，修正X坐标
X_1 = 0:rangeCell:(RFFTNum/2-1)*rangeCell;
X_Half = fliplr(X_1);       %向量翻转
V = ( -MTDFftNum/2+1 )*velocityCell:velocityCell:(MTDFftNum/2*velocityCell);    %Y轴坐标
% Z = (-FFTNumdbf/2+1)*angleCell:angleCell:(FFTNumdbf/2)*angleCell;
%Z坐标转换成m，为了观测方便
% ZmCell = maxRange/FFTNumdbf;
% Zm = (-FFTNumdbf/2+1)*ZmCell:ZmCell:(FFTNumdbf/2)*ZmCell;
%         

%采集的数据类型不同，同一驻留取出的长度不同。complex2x
sizeofData = 2;
dxstart = 1;

%记录目标轨迹变化
gfvRmax = zeros(1,FrameNum);
%定义结构体
TARGETPARA = struct(...
                    'Rr',                0.0, ...        % 目标距离，米
                    'Vr',               0.0, ...        % 目标速度，m/s
                    'Angle',            0.0, ...        % 方位角度,度
                    'Rx',               0.0, ...        %转换成直角坐标的X向距离
                    'Raz',              0.0, ...        % 目标方位向距离（角度转距离），米
                    'Amp',              0.0  ...        % 信号幅度
                    );
%记录所有帧所有点迹
TARGOUT = struct(... 
                'ValidFlag',        0.0, ...            %检测是否有效，0：没捡到目标，1：检到目标
                'TargetNum',        0.0, ...            %检测到的目标个数
                'TargetPara',       repmat(TARGETPARA, 1, TargNumMax) ...
);
gsvTargAll  = repmat(TARGOUT, 1, FrameNum);

%% 数据路径

%输入路径的方式 0：多个文件自动读入，1：选定要读的bin文件
datapathctl = 0;    


if datapathctl == 1     %只有一个bin文件
%     fname='E:\data_dca\1022_2floor\adc_data_2人_前后.bin';
        fname='E:\data_dca\1022_2floor\center\adc_data_2人_前后.bin';
%     fname='E:\data_dca\AWR1642_0925\adc_data_22.bin';
    datafileNum = 3;       %由于读取文件夹中bin文件的序号从3开始，为了方便下面的for循环，只有一个文件的时候设置序号是3
    fileCtl = 0;
else
    % 有多个bin文件，需要分别读入
%     pathName = 'E:\data_dca\1013office_2\ceil\2T_rendong';
%         pathName = 'E:\data_dca\1013office_2\adc_1p8_2T_ren';
% pathName='E:\data_dca\1022_2floor\ceil\2renxianhou';
    pathName = 'E:\data_dca\1014Parking_2\2T_miu30_128_fu2';
    
    radarPathName = strcat(pathName,'\');
    radarfile = dir(radarPathName);
    datafileNum = length(radarfile);
    fileCtl = 1;
   
end

 for dataIdx = 3:datafileNum
     % read data
    if fileCtl == 0
        fid = fopen(fname,'rb'); 
    else
        fileName = radarfile(dataIdx).name;
        fid = fopen([radarPathName,fileName],'r');
    end

    MTIFlag = 1;    %控制是否做MTI
    ThreshFlag  = 0;    %门限选择，1：慢门限，缺省：cfar
    FigNum = 0;
    ValidNum = 0;
    TargNumTemp = 0;

    colorMatrix = [0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 0 1 1 ];%黑 红 绿 蓝 洋红 青蓝

    [row_cl,col] = size(colorMatrix);

    for frameIdx = dxstart:FrameNum

        if mod(frameIdx ,100)==0
            pause(0.01);
        end

        
        if frameIdx == 120
            pause(0.01);
        end
        
        %% read data
        %读出数据 
        sdata = fread(fid,rangeSmpNum*chirpNum*...
                rxNum*txNum * sizeofData,'int16');

        % 1642+DCA1000
        %文件长度
        fileSize = size(sdata, 1);
        %复数录取，实部虚部
        lvds_data = zeros(1, fileSize/2);
        M = rangeSmpNum;
        N = chirpNum;
        %按照Data Capture.pdf上的matlab解析方法解析数据
        %LVDS使用两个通道进行数据读取
        lvds_Lanes = 2;
        Gn = lvds_Lanes * sizeofData;

        count = 1;

        for i=1:Gn:fileSize-5

           lvds_data(1,count) = sdata(i) + 1i*sdata(i+2); 

           lvds_data(1,count+1) = sdata(i+1)+1i*sdata(i+3); 

           count = count + 2;

        end
        %采集时，分时发送T1T2，分时接收，两个chirp作为一组，loops为N
        %按照 Tx1chirp1 Rx1 2 3 4 ,Tx2chirp1 Rx1 2 3 4 ,Tx1 chirp2Rx1 2 3 4 ,Tx2
        %chirp2 Rx1 2 3 4 ...
        num_chirps = txNum*N;

        lvds_data = reshape(lvds_data, M*rxNum, num_chirps);
        %注意是.'，如果只用'表示共轭转置
        lvds_data = lvds_data.';

        cdata = zeros(rxNum,num_chirps*M);

        for row = 1:rxNum

            for i = 1: num_chirps

                cdata(row,(i-1)*M+1:i*M) = lvds_data(i,(row-1)*M+1:row*M);

            end

            data_radar_Temp = reshape(cdata(row,:),M,num_chirps);   %RXn

            if txNum==1 %一根发射天线
                RxData(:,:,row) = data_radar_Temp;
            elseif txNum == 2
                %分别取出两个发射天线对应的接收数据
                Tx1Data(:,:) = data_radar_Temp(:,1:2:num_chirps);
                Tx2Data(:,:) = data_radar_Temp(:,2:2:num_chirps);
                %组成等效虚拟天线
                RxData(:,:,row) = Tx1Data;
                RxData(:,:,row+4) = Tx2Data;
            else
            end
        end

        %画出Chirp1波形
%         EchoTemp(:) = RxData(:,1,1);
%         figure(1);
%         plot(real(EchoTemp));hold on;
%         plot(imag(EchoTemp),'r'); hold off;
%         title('Chirp1波形');
%         legend('Real','Imag');

            %-------------------  计算瞬时频率  -------------------
    %         Tx1chirp1 = RxData(:,1,1);
    %         Tx2chirp1 = RxData(:,1,1);
    %         Tx1phi = phase(Tx1chirp1);
    %         Tx2phi = phase(Tx2chirp1);
    %         t_ins = 1/sampleRate*(1:SampleM);
    %         Tx1f_ins = Tx1phi./t_ins;
    %         Tx2f_ins = Tx2phi./t_ins;
    %         FigNum=FigNum+1;figure(FigNum);
    %         plot(Tx1f_ins);hold on; plot(Tx2f_ins,'r');
    %         legend('Tx1','Tx2');
    %         title('瞬时频率');
    %        

        range_win = hamming(SampleM);   %加海明窗
        %注意两级for循环与读数的for循环不同
        for ii = 1:Rx_n
             %距离维FFT
            for jj=1:ChirpN
                ffttemp = RxData(:,jj,ii).*range_win;
                FFTOut(jj,:) = fft(ffttemp,RFFTNum);
            end
            shiftfft = fftshift(FFTOut);
            %画第一个chirp的一维dB图，观测
    %         Tempfft = shiftfft(22,:);
    %         logFftOut = log2(abs(Tempfft)/max(abs(Tempfft)));
    %         FftOut = abs(Tempfft);
    %         figure(2);
    %         plot(X_All,logFftOut);
    %         titlename=['驻留号',num2str(frameIdx),'接收天线',num2str(ii),'Chirp1的一维FFT'];
    %         title(titlename)
    %         xlabel('距离（m/s）');
    %         ylabel('幅度dB')
    %         
            %FFT结果一维排列
            FFTAll = reshape(FFTOut,1,[]);

            %MTI,去除固定目标
            if MTIFlag == 1
                for jj=3:ChirpN
                    MTIOut(jj-2,:) = FFTOut(jj,:) + FFTOut(jj-2,:) -2*FFTOut(jj-1,:);
                end
                MTD_win = ChirpN - 2;
            else
                MTIOut = FFTOut;
                MTD_win = ChirpN;
            end

            vel_win = hamming(MTD_win);   %加海明窗
            %MTD，速度维FFT，用于测速
            for kk = 1:RFFTNum
                ffttemp1 = MTIOut(:,kk).*vel_win;
                MTDOutTemp(kk,:) = fft(ffttemp1,MTDFftNum);
            end
            ShiftMtd = fftshift(MTDOutTemp);

            %注意，MTDOut是M*N矩阵
            SampleMOut = RFFTNum;
            if SampleMOut == RFFTNum
                X1 = X_All;
            else
                X1 = X_Half;
            end
            MTDOut = ShiftMtd(1:SampleMOut,:);

            absMTDOut = abs(MTDOut)';
    %         figure(3);
    %         mesh(X1,Y,absMTDOut); 
    %         view(0,90);     %俯视图
    %         titlename=['驻留号',num2str(frameIdx),'接收天线',num2str(ii),'的MTD结果'];
    %         title(titlename)
    %         
    %         xlabel('距离m');
    %         ylabel('速度m/s');
    %         pause(0.01);

            absMTDOutAll(ii,:,:) = absMTDOut;
            FFT2D(ii,:,:) =MTDOut;
        end     %end of for ii = 1:Rx_n

        %通道累加，不同接收天线的相同点幅值累加  ,对累加后的幅值矩阵做cfar检测
        AmpOutAll(:,:)= sum(absMTDOutAll,1)/double(Rx_n);

        figure(4);
        mesh(X1,V,AmpOutAll); 
    %     mesh(AmpOutAll); 
        view(0,90);     %俯视图
        titlename=['驻留号',num2str(frameIdx),'的通道累加幅值矩阵'];
        title(titlename)
        xlabel('距离m');
        ylabel('速度m/s');
        pause(0.01);


            %------------------cfar 产生门限------------------
            %结构体定义
            %caCfarParaIn = struct('gap', {1,6}, 'width', {1,6}, 'cfarType', {1,6}, 'factor', {1,6}, 'offset', {1,6}, 'boundary', {1,6});
            %cfar参数
            CfarN = 2;          %保护单元个数
            CfarM = 8;         %输入滑动窗长度,积累窗口数
            CfarFactor = 4;     %门限因子
            CfarOffset = 40;   %门限直流偏置量
            %两个一维相关区幅度最大值对应位置相差不超过Rdet个距离门时可以进行二维相关
            Rdet = 5;
            %认为是一个目标的速度相关区长度
            %Vdet = 4;

            %结构体参数
            caCfarParaIn.gap 		= CfarN;        %保护单元个数
            caCfarParaIn.width 		= CfarM;        %输入滑动窗长度
            caCfarParaIn.factor 	= CfarFactor;	%门限因子
            caCfarParaIn.offset		= CfarOffset;	%门限直流偏置量
            caCfarParaIn.cfarType	= 0;            %输入滑动平均处理方式0:滑动平均;1:滑动平均选大; 2:滑动平均选小,3:慢门限
            caCfarParaIn.boundary	= 1;            %边界处理方式0:固定值;1:单边门限值

            %所有chirps每个距离单元生成cfar门限
            Vlg = zeros(MTDFftNum,SampleMOut);
            for kk1 = 1:MTDFftNum
                %生成缓冲
                buffData=zeros(1,SampleMOut);
                if ThreshFlag == 1
                    %慢门限
                    %设慢门限计算门限使用的距离单元个数为
                    ThreshCountNum = SampleMOut/4;
                    SlowCoeff = 120;
                    [Threshold(kk1,:),reFlag(kk1)] = SlowThresh(AmpOutAll(kk1,:), SampleMOut, SlowCoeff, ThreshCountNum);
                else
                    %计算cfar门限值
                    [Threshold(kk1,:),reFlag(kk1)] = caCfar(caCfarParaIn, AmpOutAll(kk1,:), buffData, SampleMOut);
                end
                %大于门限填1
                Vlg(kk1,:) = AmpOutAll(kk1,:) > Threshold(kk1,:);

            end
            aa = find(Vlg > 1e-3); %find查找的时候，先列后行
            %------------------相关
            % 一维相关 距离维限制0-5
            [TotalTgtNum,Tinfo2]=D1correlation2(Vlg,Threshold,caCfarParaIn.factor,AmpOutAll,MTDFftNum,SampleMOut,Rdet);

            % 二维相关 速度维限制0-5
            %struct Plot
            %{
            %	V;VStart;VEnd;
            %	Sigmags;RMax;VMax;
            %	R;RStart;REnd;
            %	Azimuth;Time;Flag
            %}
            [Plot,CurrentTarNum]=D2correlation2_dsp(Tinfo2,MTDFftNum,TargNumMax,MTDFftNum);

            %取出二维相关后各目标距离
            if CurrentTarNum > 0
                %取出的点迹信息按照最大值方法还是中心法，0:重心，1:最大值
                Calflag = 1;
               [fvRIndex ,fvVIndex, AmpAll] = GetAllTargCoordinate(Calflag , Plot );

            end
            fvAmax = zeros(SampleMOut,MTDFftNum);

            for ii = 1:CurrentTarNum
                fvAmax(fvRIndex(ii),fvVIndex(ii)) = AmpAll(ii);
            end


    %         if CurrentTarNum > 0
    %             figure(8);
    %             mesh(fvAmax');
    %             title('Cfar检测结果（相关后）')
    %             xlabel('距离维（采样点）');
    %             ylabel('多普勒维（采样点）');
    %         end
    %         
            %保存各个接收检测到的目标个数
            ivTargNum = CurrentTarNum;
            %-----------------搜索目标结束

            pause(0.01);

        %在FFT2D中取出各个目标的复数值
        for ii = 1:ivTargNum

            %分别得到距离维速度维坐标
            RIndex = fvRIndex(ii);
            VIndex = fvVIndex(ii);
            %取出各接收对应点信息，用于后续不同接收天线之间FFT求角度
            for jj = 1:Rx_n
                cvTarg(jj,ii) = FFT2D(jj,RIndex,VIndex);
            end  

        end

        %不同通道的相同目标，作FFT，等效DBF
        fvDBFAmax = zeros(1,SampleMOut);
        %数组定义，初始化
        fvTargTheta = zeros(1,FFTNumdbf);
        fvTargR = zeros(1,TargNumMax); 
        Ymeter = zeros(1,TargNumMax); 
        TargInfo = zeros(SampleMOut,FFTNumdbf);

        TargNumOut = 0;
        for ii = 1:ivTargNum
            DBFOut = fftshift(fft(cvTarg(:,ii),FFTNumdbf));
    %         DBFOut = fft(cvTarg(:,ii),FFTNumdbf);

            [Amax,iChnNum] = max(abs(DBFOut));
            Rtemp= X1(fvRIndex(ii));

            %距离为负，无效，去除，只保留正距离的目标
            if  Rtemp >= 0 
                %输出目标个数
                TargNumOut = TargNumOut +1;
                %目标方位角
                fw = (iChnNum-FFTNumdbf/2-1)/FFTNumdbf;
                fvTargTheta(TargNumOut) = asin(fw * lamda / dPerRx) * 180 / pi;
            
                %目标距离
                fvTargR(TargNumOut) = Rtemp;
                %转化为直角坐标
                Xmeter(TargNumOut) = cos(fvTargTheta(TargNumOut)*pi/180)*fvTargR(TargNumOut);
                %方位方向的距离
                Ymeter(TargNumOut) = sin(fvTargTheta(TargNumOut)*pi/180)*fvTargR(TargNumOut);
                %目标速度
                fvTargV(TargNumOut) = V(fvVIndex(ii));
                %幅度
                fvDBFAmax(TargNumOut) = Amax;
                TargInfo(fvRIndex(ii),iChnNum) = Amax;
                %所有目标点迹保存
                gsvTargAll(frameIdx).TargetPara(TargNumOut).Rr = X1(fvRIndex(ii));
                gsvTargAll(frameIdx).TargetPara(TargNumOut).Vr = V(fvVIndex(TargNumOut));
                gsvTargAll(frameIdx).TargetPara(TargNumOut).Angle = fvTargTheta(TargNumOut);
                gsvTargAll(frameIdx).TargetPara(TargNumOut).Rx = Xmeter(TargNumOut);
                gsvTargAll(frameIdx).TargetPara(TargNumOut).Raz = Ymeter(TargNumOut);
                gsvTargAll(frameIdx).TargetPara(TargNumOut).Amp = Amax;
            end

        end
        %回送标志和目标个数
         if TargNumOut>0
            gsvTargAll(frameIdx).ValidFlag = 1;
            gsvTargAll(frameIdx).TargetNum = TargNumOut;
         else
             gsvTargAll(frameIdx).ValidFlag = 0;
            gsvTargAll(frameIdx).TargetNum = 0;
         end

        %最大目标
        [Ampmax,maxindex] = max(fvDBFAmax);
        gfvRmax(frameIdx) = fvTargR(maxindex);


        if SampleMOut == RFFTNum
            %负距离无效
            TargInfoOut = TargInfo(SampleMOut/2+1:SampleMOut,:);
    %         TargInfoOut = TargInfo(1:SampleMOut,:);
        else
            TargInfoOut = TargInfo;
        end
    %     figure(6);
    % %     subplot(2,1,1);
    %     handle = mesh(Z,X_1,TargInfoOut);
    % %     handle = mesh(TargInfoOut);
    %     view(0,90);
    %     title('目标运动信息');
    % %     text1 = text(fvTargTheta(1:TargNumOut),fvTargR(1:TargNumOut)+0.4,num2cell(fvTargV(1:TargNumOut)),'Color','red');
    %     xlabel('方位角度（度）')
    %     ylabel('目标距离（米）')
        if TargNumOut>0
                %有效个数，用于画图
                ValidNum = ValidNum + 1;
        end

        %保存当前周期检测结果
        if TargNumOut>0
            TargNumTemp = TargNumOut;
            YOutTemp = Ymeter(1:TargNumOut);
            XOutTemp = Xmeter(1:TargNumOut);
        end
        if mod(frameIdx,1)==0
            if TargNumTemp > 0
                figure(7);
                %画图颜色
                Cindex = floor(ValidNum/1);
                colorIndex = mod(Cindex,row_cl) + 1;

                h(ValidNum) = plot(YOutTemp(1:TargNumTemp),XOutTemp(1:TargNumTemp),'*','color',colorMatrix(colorIndex,:));
%                 h(ValidNum) = plot(YOutTemp(1:TargNumTemp),XOutTemp(1:TargNumTemp),'*b');
%                  text1 = text(YOutTemp(1:TargNumTemp),XOutTemp(1:TargNumTemp)+0.04,num2cell(fvTargV(1:TargNumTemp)),'Color','red');
                
                %限定坐标
%                 axis([-5 5 0 10]);
                title('目标信息');
                xlabel('方位向距离（米）')
                ylabel('目标距离（米）')
                hold on;
                pause(0.01);
               %随时间，旧图删除，也可以用delete命令
        %         if ValidNum > 40
        %             set(h(ValidNum - 40),'visible','off'); %不显示
        %         end
        %          
            end
            %画图后，清空临时缓冲
            TargNumTemp = 0;
        end
    %     mesh(Zm,X_1,TargInfoOut);
    %     view(0,90);
    %     title('目标运动信息');
    %     xlabel('方位向距离（米）')
    %     ylabel('目标距离（米）')
    %     waitforbuttonpress;


    %20201010临时添加，两维距离差都在0.5m以内的的点认为是同一个目标
        %对点迹进行合并，相同目标回送幅度较大的点
    %    RZ0P5To1Targ;


        pause(0.01);

    end
    
 end
fclose(fid);
%保存所有帧点迹，用于运动轨迹跟踪处理
% save 'TargInfoAll' gsvTargAll




