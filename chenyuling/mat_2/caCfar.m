function[y, reflag] = caCfar(caCfarparaIn, x, buff, n) 

%初始化计算变量
k1=0;
reFlag=0;
width=0;gap=0;cfarType=0;offset=0;
factor=0;
y(1) = 0;

%滑动平均函数参数结构变量赋给内部变量
width=caCfarparaIn.width;
gap=caCfarparaIn.gap;
cfarType=caCfarparaIn.cfarType;
factor=caCfarparaIn.factor;
offset=caCfarparaIn.offset;
boundary=caCfarparaIn.boundary;
%如果处理数据长度不满足要求
    
if (n<=(2*width+2*gap+1) || n<=1) || width<=0
    %异常处理,置函数返回值为-1
    reflag=-1;
else
    %进行滑动平均处理
    for k1=1:n-width+1
        buff(k1) = (factor/width) * sum(x(k1:(k1+width-1))) + offset;
    end  
    %滑动平均处理方式选择
    switch (cfarType)
        %滑动平均
        case 0          
            %左门限为k,右门限为k+_win+2*_gap+1,对应点为k+_win+_gap
            %求和计算
            for k1=1:n-2*width-2*gap
                y(width+gap+k1) = buff(k1) + buff(k1+width+2*gap+1);
            end
            %求和结果求平均
            y = y * 0.5; 
            %标志位正常
            reflag=0;
        %滑动平均选小
        case 1          
            %左门限为k,右门限为k+_win+2*_gap+1,对应点为k+_win+_gap
            %取小计算
            for k1=1:n-2*width-2*gap
                y(width+gap+k1) = min(buff(k1),buff(k1+width+2*gap+1));
            end
            %标志位正常
            reflag=0;
        %滑动平均选大
        case 2          
            %左门限为k,右门限为k+_win+2*_gap+1,对应点为k+_win+_gap
            %取小计算
            for k1=1:n-2*width-2*gap
                y(width+gap+k1) = max(buff(k1),buff(k1+width+2*gap+1));
            end
            %标志位正常
            reflag=0; 
        %不在处理方式内的默认处理为滑动平均处理,函数返回值为-2
        otherwise            
            %左门限为k,右门限为k+_win+2*_gap+1,对应点为k+_win+_gap
            %求和计算
            for k1=1:n-2*width-2*gap
                y(width+gap+k1) = buff(k1) + buff(k1+width+2*gap+1);
            end
            %求和结果求平均
            y = y * 0.5;            
            %标志位正常
            reflag=-2; 
    end %end switch (cfarType)
    
    %左右边界处理
    if boundary %左右边界按照单边处理
        for k1=1:width+gap
            y(k1) = buff(1+gap+k1);
            y(n-k1+1) = buff(n-width-gap-k1+1);  
        end
    else %左右边界按照最后一个有效固定值处理
        for k1=1:width+gap
            y(k1) = y(width+gap+1);
            y(n-k1+1) = y(n-width-gap);  
        end       
    end
end


