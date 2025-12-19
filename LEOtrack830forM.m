%% LEO 跟踪代码  许哲锴
%2020 2 10  加入各阶数的PLL,FLL，修改注释
%2025.5.12  加入LKF EKF UKF
%2025 5 20  加入不需要剥离电文的EKF(Rkf还需要推导)暂时不能用 mode4
%2025.6.7   加入HDC-FLL
%2025 6 10  加入CSK解调 mode5 不做任何调整
%2025 7.5   加入计算可视范围内的每1ms多普勒 修改多普勒更新的模型
%2025.7.8   加入LKF EKF 的外推算法，并对外推条件进行限制,mode8 和mode 9
%2025.7.11  加入LDPC的编码，修改高速电文，4个符号帧头+6个任意符号 循环播发
%2025 7.13  加入LDPC解码函数LDPC_decode_main(LLR) 存在一定问题 该函数不用了
%2025 7 14  重构高速电文
%2025 7 15  新增LDPC解码函数 统计BER 跟踪误差函数
%2025 7 15  mode5保持相干时本地与接受相同
%2025 7 16  修改部分矩阵的拼接，优化速度
%2025 7 16  新增解调判决条件
%2025 7 17  修复载波公式的BUG，另外将跟踪误差，解调数据保存至Excel中
%2025 8 30  测试跟踪概率
clear all
clc;
%% Monte-Carlo harness
num_mc = 5;                                  % 运行次数
seeds  = 1:num_mc;                            % 或者 randi(1e9,1,num_mc)
doppler_fail_threshold = 250;                % 失败阈值（Hz），可改
consecutive_ms_for_fail = 100;   % 连续 N ms 超阈才算失锁; 若按论文0.5s, 这里=500
% 汇总表预分配（12列：seed/CN0/mode/fail_ms/...）
mc_log = array2table(nan(num_mc,12), 'VariableNames', ...
 {'seed','CN0','mode','fail_ms','BERcoh','BERnoncoh','BERldpc', ...
  'RMSEphi','RMSEdop','CN0mean','runtime_s','success'});

for mc = 1:num_mc
    fprintf('==== MC %d/%d ====\n', mc, num_mc);
    rng(seeds(mc));                 % 每次不同随机种子
    %% 计算可视范围内的每1ms多普勒
    % 一些常数和参数设置
    vs = 7618.81;     % Satellite velocity in m/s
    re = 6368000;     % Earth's radius in meters
    rs = 6868000;     % Satellite radius in meters
    w = 0.00110932;
    G = 6.67430e-11;  % Gravitational constant (m^3 kg^-1 s^-2)
    M = 5.972e24;     % Earth's mass (kg)
    
    % 角度范围(并不是仰角，指的是卫星与地心的夹角)
    anger=80;
    x = linspace(anger*pi/180, pi - anger*pi/180, 10000);  % Angle range for satellite motion  80
    r_sine=rs*cos(anger*pi/180)/sqrt(re^2+rs^2-2*re*rs*sin(anger*pi/180));
    elevation_rad=acos(r_sine);
    elevation_degree=elevation_rad/(2*pi)*360; %仰角 degrees
    
    % 计算多普勒频偏 yv 和多普勒一阶变化率 ya
    yv = (1530 / 300) * vs * re .* cos(x) ./ sqrt(re^2 + rs^2 - 2 * re * rs * sin(x));
    ya = -(1530 / 300) * w .* vs .* re .* (re .* rs .* sin(x) .* sin(x) - (re .* re + rs .* rs) .* sin(x) + re .* rs) ./ (re .* re + rs .* rs - 2 * re .* rs .* sin(x)).^(3/2);
    
    % 计算卫星轨道周期 (运动时间)
    a = rs;  % 半长轴
    T = 2 * pi * sqrt(a^3 / (G * M));  % 轨道周期 (秒)
    T_minutes = T / 60;  % 转换为分钟
    
    % 计算时间与角度的关系
    time = (x / (2 * pi)) * T;  % 对应的时间
    escpae=time(end)-time(1);%卫星的过境时长 单位 s
    % 为了每 1ms 计算一次多普勒，我们用 1ms 的时间步长来离散化时间
    dt = 1e-3;  % 1ms = 0.001秒
    time_steps = time(1):dt:(time(end)+10);  % 从0到轨道周期，步长为1ms  看time的头和末尾
    num_steps = length(time_steps);  % 时间步数
    
    % 根据每1ms计算时间步长，重新计算多普勒频率和变化率
    yv_new = zeros(1, num_steps);
    % ya_new = zeros(1, num_steps);
    
    for i = 1:num_steps
        % 在每个时间点，计算对应的角度
        angle = 2 * pi * (time_steps(i) / T);  % 时间对应的角度
        % 计算对应的多普勒频偏和变化率
        yv_new(i) = (1530 / 300) * vs * re * cos(angle) / sqrt(re^2 + rs^2 - 2 * re * rs * sin(angle));
        %     ya_new(i) = -(1530 / 300) * w * vs * re * (re * rs * sin(angle)^2 - (re^2 + rs^2) * sin(angle) + re * rs) / (re^2 + rs^2 - 2 * re * rs * sin(angle))^(3/2);
    end
    
    %配置参数
    runtime = round(escpae)+1;                                             %程序运行时间，单位秒
    filename = 'LKF.txt';                                       %打印到文件
    format long
    %%%%%%% 定义参数部分：
    speedoflight = 299792458;
    Freq_sample = 16.369e6;                                     %%%%%%%%采样频率
    Fc = 4092000;                                              %%%%%%%% 低中频无多普勒频偏载波频率
    Freq_code = 2.046e6;                                       %%%%%%%% 扩频码频率
    % acce = [200  0 0 0];                                         %%%%%%%% 径向加速度
    % jerk=[-10 0 0 0];
    Freq_IF = [Fc+yv_new(1)  Fc+1000   Fc+2500  Fc-1000];           %%%%%%%% 实际低中频载波频率
    Frf = 1522224000;                                          %%%%%%%% 射频频率
    codedopcoef = Frf/Freq_code;                               %射频频率与扩频码频率之比
    codelength = 2046;                                         %码长
    bitperiod = 4;                                            %bit周期，单位ms
    % freq_dot = acce/speedoflight*Frf;                                   %%%%%%%% 根据公式 f = v / c * Frf得频率变化率 f' = acce / c * Frf
    % freq_code_dot = acce/speedoflight*Freq_code;                        %%%%%%%% 同理可得伪码频率变化率为 f_code' = acce / c * Freq_code
    % freq_dot_2rd=jerk/speedoflight*Frf;                                  %二阶变化率
    % freq_code_dot_2rd=jerk/speedoflight*Freq_code;
    chipcounter = [0 0 0 0];                                   %%%%%%%% 码片计数器，四颗卫星
    carrier_dco_delta = (Freq_IF)/Freq_sample;                %%%%%%%% 载波以cos(2*pi*fc*t)的形式在信号中，离散化后变成cos(2*pi*fc*n*Ts),对于零中频信号则为(Freq_IF - Fc) / Freq_sample
    code_dco_delta = ((Freq_IF-Fc)/codedopcoef+Freq_code)/Freq_sample;%%%%%%%% 同上
    carrier_dco_phase = [0 0 0 0];                             %%%%%%%% 初始化载波相位，四颗卫星
    code_dco_phase=[0 0 0 0];                                  %%%%%%%% 初始化伪码相位，四颗卫星
    mscount = zeros(1,4);                                      %%%%%%%% 毫秒计数器，四颗卫星
    
    bitnumlow = zeros(1,4);                                     %%%%%%%% 比特计数器，四颗卫星(低速：250bps)
    bitnumhigh = zeros(1,4);                                    %%%%%%%% 比特计数器，四颗卫星(高速：3000bps)
    gt = zeros(1,4);                                             %高速低速选通器,0低速选通，1高速选通
    
    
    %%%%%%%% 调用码发生器函数生成四颗卫星的伪码
    codelow=zeros(4,2046);
    codehigh_initial=zeros(4,2046);
    %%调用码发生器函数，两个，低速高速
    for i=1:1:4
        [codetemphigh,codetemplow]=CAgenerate(i);
        codehigh_initial(i,:)=codetemphigh;
        codelow(i,:)=codetemplow;
    end
    
    %%%%%%%% 量化判决前通过一个低通滤波器，定义如下：
    filterorder = 21;                                                  %%%%%%%% 阶数为20，此处表示有21个点
    result_input = zeros(1,filterorder);                               %%%%%%%% 通过滤波器的输入数据数组，长度等于filterorder
    %%% 生成低通滤波器：选用Kaiser窗FIR滤波器，Beta值为0.5，阶数为20，采样率30.69MHz，截止频率1.1MHz，系数如下：
    % lowpassfilter = [0.021296855511562 0.0276752982618497 0.0341164198702814 0.0404241078806896 0.0463987613744401...
    %                 0.051845537789563 0.0565825791873747 0.060448817384489 0.0633109732940896 0.0650694024346712...
    %                 0.0656624940219796 0.0650694024346712 0.0633109732940896 0.060448817384489 0.0565825791873747...
    %                 0.051845537789563 0.0463987613744401 0.0404241078806896 0.0341164198702814 0.0276752982618497...
    %                 0.021296855511562
    %                 ];                               %%%%%%%% 1*21的滤波器系数
    %%% 生成低通滤波器：选用Kaiser窗FIR滤波器，Beta值为0.5，阶数为20，采样率16.368MHz，截止频率1.1MHz，系数如下：
    lowpassfilter = [  -0.024878413973009  -0.019403240492844  -0.008442779911688   0.007679140518129...
        0.027927884076498   0.050642413108743   0.073714613080674   0.094841094821135...
        0.111813064702503   0.122803007209291   0.126606433721137   0.122803007209291...
        0.111813064702503   0.094841094821135   0.073714613080674   0.050642413108743...
        0.027927884076498   0.007679140518129  -0.008442779911688  -0.019403240492844...
        -0.024878413973009
        ];                               %%%%%%%% 1*21的滤波器系数
    
    %%% 将滤波器系数归一化，使得数据通过滤波器后幅值不变
    H = sum(lowpassfilter.*lowpassfilter);
    lowpassfilter = lowpassfilter/sqrt(H);
    %% 电文生成
    codehigh = codehigh_initial;%高速码设两个，为移位做准备
    xk = zeros(1,4);   %6bit拼接符号
    %% low bit 250bps
    % 帧头（24比特）
    FrameSynFlag = [1 1 1 0 0 0 1 0 0 1 0 0 1 1 0 1 1 1 1 0 1 0 0 0];
    
    % 尾部（66比特），重复11次得到726比特
    Tail1 = [0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,...
        0,0,1,1,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,0];
    Tail = repmat(Tail1, 1, 11);  % 726比特
    
    % 构建完整帧（750比特）
    Frame = [FrameSynFlag Tail];
    
    % 生成 500*1 个帧，共 375000*1 比特
    BITlow = repmat(Frame, 1, 500*1);
    BITlow = BITlow*2-1;
    
    
    %高速电文 3000bps
    load('LDPC200_100high');  % 加载 LDPC 编码比特 resultbin1，长度应为 1200
    
    BIThigh = repmat(resultbin1, 1, 1000);  % 1200 × 1000 = 900000 bits
    
    %%%%%%%% 设定判决门限：
    sigma = 1;
    th = 1;
    BITHIGHTRANS=BIThightrans(BIThigh);
    %%%%%%%% 根据公式：(1)sigma^2 =  *N0/2; (2)sigma0 = N0*B;
    %%%%%%%% (3)CNR - 10*lgB = 10*lg(A^2/2/sigma0^2) 得
    %%%%%%%% A = sigma * sqrt(2 * 10^(CN0/10)/fs)   ********
    CN0 = 34;
    A(1) = sigma* sqrt(2*10^(CN0/10)/Freq_sample);%
    A(2) = sigma* sqrt(2*10^(CN0/10)/Freq_sample);%
    A(3) = sigma* sqrt(2*10^(CN0/10)/Freq_sample);%
    A(4) = sigma* sqrt(2*10^(CN0/10)/Freq_sample);%
    
    n = 0;
    flag = 0;
    num = 0;
    samplenum = 0;
    sec = 0;
    
    
    %跟踪参数配置以及累加数组初始化
    mode = 8;                                          %设定模式 . ==0 pll   ==1 fll   ==2 LKF   ==3 EKF剥离电文   ==4 EKF不剥离电文  ==5 不作任何调整，用于CSK解调  ==6 UKF剥离电文
    %mode==8  LKF 1MS预测2ms矫正
    %mode==9  EKF 1MS预测2ms矫正剥离电文
    
    AcqSatNUM = 1;                                     %设定要跟踪的卫星PRN编号
    LocalDop = Freq_IF(AcqSatNUM)-Fc+20;                %初始载波多普勒，注意要跟设定值比较接近
    %mode==5时，重新计算之间的步进以覆盖，确保多普勒为0。
    if(mode==5)
        yv_new = zeros(1, num_steps);
        LocalDop=0;
        carrier_dco_delta(AcqSatNUM) = (Fc)/Freq_sample;
        code_dco_delta(AcqSatNUM) = ((Fc-Fc)/codedopcoef+Freq_code)/Freq_sample;
    end
    
    [codetemphigh,codetemplow]=CAgenerate(AcqSatNUM);
    LocalCode = codetemplow;                                       %初始跟踪通道扩频码
    for cskn=0:1:63
        Localhightraverse(cskn+1,:)=circshift(codetemphigh,cskn);
    end
    demodulation=1;%解调开关
    if(CN0<35)
        demodulation=0;
    end
    if(CN0>30)
        jugujitime=80;
    else
        jugujitime=1000;
    end
    fprintf('当前选择跟踪模式=%d   设定载噪比= %f    运行时长=%d 秒  仰角=%.4f ° 解调模式开启1，关闭0 ：目前为%d\n', mode,CN0,runtime,elevation_degree,demodulation);
    
    
    code_dco_phaseL = 0;                               %本地复现码相位
    carrier_dco_phaseL = 0;                            %本地复现码相位
    % code_dco_phaseL2 = 0 ;%本地复现码相位，不归一化到0-1，后面没有用到该量
    code_dco_deltaL = ((LocalDop)/codedopcoef+Freq_code)/Freq_sample;     %本地码频率步进
    carrier_dco_deltaL = (Fc + LocalDop)/Freq_sample;                     %本地载波频率步进
    chipcntL = 0;                                     %本地码片计数
    mscountL = 0;                                     %本地ms计数
    bitnumL = 0;                                      %本地bit计数
    mstime=0;                                         %用于统计跟踪误差
    LocalCorreArray = [-1.5:0.1:1.5];                 %本地码搜索码片范围和码片间距
    LengthCorreArray = length(LocalCorreArray);       %相关器数目
    mid = (LengthCorreArray+1)/2;                     %P路索引
    corrsum = zeros(1,LengthCorreArray);              %所有支路累加值
    corrcohsum_carr = zeros(1,LengthCorreArray);      %相干载波环累加值
    corrcohsum_code = zeros(1,LengthCorreArray);      %相干码环累加值
    corrnoncohsum_code = zeros(1,LengthCorreArray);   %非相干码环累加值
    % corrcohsum_carr_FLL=zeros(1,LengthCorreArray);    %FLL累加值
    corrcohsum_carr_FLL_first=0;                               %存放这一次FLL的累加值
    corrcohsum_carr_FLL_second=0;                               %存放这上次FLL的累加值
    dump = 0;                                         %为1时，表示一个周期的累加结束，一般情况一个周期表示1个扩频码周期，B1C等例外
    localn = 0;                                       %本地采样点计数
    FLL_bitcount=0;                                        %FLL用于记录当前相干到第几个比特周期
    %%
    %码环路参数配置
    %%%%%%%%%%%本地码环参数%%%%%%%%%%%%%
    T = 0.001;                                        %单位为秒，表示1ms更新一次
    cd = 0;                                           %码鉴相器值
    cd_1 = 0;                                         %码鉴相器上一个时刻的值
    vo = 0;                                           %压控振荡器输出
    vo_1 = 0;                                         %压控振荡器上一个时刻的值
    codeloop_cohtime = 4;                             %码环相干时间，2表示2个累加周期，一般是2ms，要具体分析
    codeloop_noncohtime = 2;                         %码环非相干累加次数
    codeloop_noncohcnt = 0;                           %码环非相干累加计数器
    bn_code = 1;                                      %码环带宽1
    cd1 = 8/3*bn_code;                                %二阶码环参数一
    cd2 = cd1*cd1/2*T*codeloop_cohtime*codeloop_noncohtime/2;   %二阶码环参数二，K2
    %%
    %载波环路PLL参数配置
    %%%%%%%%%%%%%%%%%%%%%%本地载波环参数%%%%%%%%%%%%%%%%%%%%%%%
    carrierlooporder = 3;                             %载波环阶数配置，2表示2阶环路，3表示3阶环路
    pd = 0;                                           %载波鉴相器输出
    pd_1 = 0;                                         %上一次鉴相器输出
    pd_2 = 0;                                         %再上一次的鉴相器输出
    pd_3=0;
    vn = 0;                                           %压控振荡器输出
    vn_1 = 0;                                         %上一次压控振荡器输出
    vn_2 = 0;                                         %再上一次压控振荡器输出
    vn_3=0;
    carrierloop_cohtime = 4;                          %环路相干积分时间，单位ms
    bn_carr = 25;                                     %载波环路带宽
    %二阶
    cp1_2nd = 8/3*bn_carr;                            %二阶环路系数1
    cp2_2nd = cp1_2nd*cp1_2nd/2;                      %二阶环路系数2
    %三阶
    cp1_3rd = 60/23*bn_carr;                          %三阶环路系数1
    cp2_3rd = cp1_3rd*cp1_3rd*4/9;                    %三阶环路系数2
    cp3_3rd = (cp1_3rd*cp1_3rd*cp1_3rd)*2.0/27.0;     %三阶环路系数3
    %四阶
    cp1_4th=64/27*bn_carr;
    cp2_4th=(cp1_4th^2/2);
    cp3_4th=(cp1_4th^3/8);
    cp4_4th=(cp1_4th^4/64);
    
    %%
    %载波环路FLL参数配置
    %%%%%%%%%%%%%%%%%%%%%%本地锁频环参数%%%%%%%%%%%%%%%%%%%%%%%
    carrierloopFLLorder = 2;                          %锁频环阶数配置，2表示2阶环路，3表示3阶环路
    fd = 0;                                           %载波鉴相器输出
    fd_1 = 0;                                         %上一次鉴相器输出
    fd_2 = 0;                                         %再上一次的鉴相器输出
    fd_3=0;
    fo_n = 0;                                           %压控振荡器输出
    fo_n_1 = 0;                                         %上一次压控振荡器输出
    fo_n_2 = 0;                                         %再上一次压控振荡器输出
    fo_n_3 = 0;                                         %再上一次频率输出
    carrierloop_fll_cohtime = 4;                       %FLL环路的相干积分时间取1个比特周期，单位ms
    carrierloop_fll_real_cohtime=carrierloop_fll_cohtime/2; %FLL单次相干积分时长 需要根据这个参数值设定初始载波多普勒偏差！！！
    %1个电文比特的周期分成前一半和后一半，相当于单次积分时间为carrierloop_fll_cohtime/2
    %atan的鉴别范围为±1/(4*carrierloop_fll_cohtime/2) Hz
    %即为±1/(4*carrierloop_fll_real_cohtime)Hz
    carrierloop_fll_noncohcnt = 0;                       %当前用了几个周期的电文进行相干
    carrierloop_fll_noncohtime =2;                      %总计使用几个电文周期进行相干
    bn_carr_fll =3;                                     %载波环路带宽 38dbHZ时bn10可以跟踪 5
    cross_dot=0;                                         %叉积点击初始化
    cf1_1st = 4.0 * bn_carr_fll * T * carrierloop_fll_cohtime*carrierloop_fll_noncohtime/2.0;
    %2阶FLL系数
    c1 = 8.0/3.0 * bn_carr_fll;
    c2 = (c1*c1)/2.0;
    cf1_2nd = c1 * (T*carrierloop_fll_cohtime*carrierloop_fll_noncohtime)/2.0;
    cf2_2nd = c2 * ((T*carrierloop_fll_cohtime*carrierloop_fll_noncohtime)/2.0)^2;
    %3阶FLL系数
    c1 = 60.0/23.0 * bn_carr_fll;
    c2 = (c1*c1)*4.0/9.0;
    c3 = (c1*c1*c1)*2.0/27.0;
    cf1_3rd = c1 * (T*carrierloop_fll_cohtime*carrierloop_fll_noncohtime)/2.0;
    cf2_3rd = c2 * ((T*carrierloop_fll_cohtime*carrierloop_fll_noncohtime)/2.0)^2;
    cf3_3rd = c3 * ((T*carrierloop_fll_cohtime*carrierloop_fll_noncohtime)/2.0)^3;
    
    %%
    %跟踪数组初始化
    lenthofmatrix=(runtime+10)*1000;%定义矩阵长度
    FCODE=zeros(1,lenthofmatrix);
    FCARR=zeros(1,lenthofmatrix);
    
    Earr = zeros(1,lenthofmatrix); %统计E路
    Parr=  zeros(1,lenthofmatrix); %统计P路
    Larr=  zeros(1,lenthofmatrix); %统计L路
    
    loopcnt = 0;               %环路计数
    cntL = 0;                  %未使用
    
    PD=    zeros(1,lenthofmatrix);
    PDFLL= zeros(1,lenthofmatrix);
    CN0L=  zeros(1,lenthofmatrix);
    codeerr=zeros(1,lenthofmatrix);
    carrerr=zeros(1,lenthofmatrix);
    DOPERRORre=zeros(1,lenthofmatrix);
    
    flag = 1;                  %为1表示使用相位延展
    testopen = 1;
    M2=0;%用于在载噪比估计
    M4=0;
    corrsumforcskdemodulate=zeros(1,64); %用于CSK解调
    
    highbitcsknoco_symbol=zeros(1,lenthofmatrix);
    highbitcskco_symbol=zeros(1,lenthofmatrix);
    %LDPC相关运算数组
    max_decode_times = 5000;  % 预估最大解码次数（如：300秒×10次/秒 = 3000次）
    decode_matrix = zeros(200, max_decode_times);  % 每列存一次 decode
    decode_index = 1;
    highbitcskcoforLDPC=[];
    highbitcskcosysafteerLDPC=[];
    
    tt1=0;
    acccom=0;
    n1=0;
    avgtime=1000;
    disp("已完成初始化，开始跟踪");
    %输出提示信息：
    if(mode==0)
        fprintf('当前选择跟踪模式为PLL，环路阶数 % d, 带宽 % d   更新时长 %d  ms,有效积分时长 % d ms\n', ...
            carrierlooporder,bn_carr,carrierloop_cohtime,carrierloop_cohtime/2);
    elseif(mode==1)
        fprintf('当前选择跟踪模式为HDC-FLL，环路阶数 % d, 带宽 % d   电文组数 % d ,更新时长%d ms，单次相干积分时长 % d ms，实际单次有效积分时长 % d ms\n', ...
            carrierloopFLLorder,bn_carr_fll,carrierloop_fll_noncohtime,carrierloop_fll_noncohtime*4,...
            carrierloop_fll_real_cohtime,carrierloop_fll_real_cohtime/2);
    end
    
    t0=tic;  % 开始计时
    timeold=0;
    t=0;
    h = waitbar(0, '跟踪处理中...', 'Name', '跟踪进度');
    success=1;
    waitbar(0 / runtime, h, sprintf('第 %d/%d 秒数据，已耗时 %.1f 秒', 0, runtime, 0));
    bad_ms = 0;  % 连续超阈(ms)计数器，用于“连续N毫秒超阈才失锁”的口径
    %% 信号生成及跟踪
    while (1)
        %%
        %该部分为信号跟踪PLL+FLL+码环路鉴相，调整模块
        if(dump == 1)
            loopcnt = loopcnt + 1;
            dump = 0;
            if(mod(loopcnt,2)==0)
                cntL = cntL + 1;
            end
            
            if(testopen)
                Earr(loopcnt) = corrsum(mid-5);
                Larr(loopcnt) = corrsum(mid+5);
                Parr(loopcnt) = corrsum(mid+0);
            end
            if(demodulation)
                %统计CSK解调数学 判断条件是 mstime/2==0?
                if(mod(mstime,2)==0)
                    %非相干解调
                    [~, maxIndex] = max(abs(corrsumforcskdemodulate));
                    
                    highbitcsknoco_symbol(cntL)=(maxIndex-1);
                    %相干解调
                    [~, maxIndex] = max(abs(real(corrsumforcskdemodulate)));%对real取绝对值，防止载波相位存在180°的模糊
                    
                    highbitcskco_symbol(cntL)=(maxIndex-1);
                    
                    % 暂存每一个符号的64组相干值
                    highbitcskcoforLDPC=[highbitcskcoforLDPC corrsumforcskdemodulate'];
                    if(mod(mstime,400)==0)
                        iter_max=30;
                        GFq=64;
                        nm=8;
                        [iter_res, decode,GF_cnt,REAL_cnt] = EMSdecoderforcsk(highbitcskcoforLDPC,iter_max, GFq, nm);
                        
                        % 存入预分配矩阵
                        decode_matrix(:, decode_index) = decode(:);  % 强制转为列向量
                        decode_index = decode_index + 1;
                        
                        highbitcskcoforLDPC=[];
                    end
                    corrsumforcskdemodulate=zeros(1,64);
                end
                %将每一个符号的64组相关值保存，当保存了200个符号时，可以进行一次纠错编码
            end
            
            corrcohsum_carr = corrcohsum_carr + corrsum;         %相干累加为了载波更新
            corrcohsum_code = corrcohsum_code + corrsum;         %相干累加为了码环更新
            if(mscountL == 0)                                    %
                mscountL1_bitperiod = bitperiod;
            else
                mscountL1_bitperiod = mscountL;
            end
            %         accum(mscountL1_bitperiod) = corrsum(mid+0);           %accum用来计算载噪比和PLD  扎带款单
            M4=M4+(real(corrsum(mid+0))^2+imag(corrsum(mid+0))^2)^2;
            M2=M2+real(corrsum(mid+0))^2+imag(corrsum(mid+0))^2;
            
            %使用FLL时，单次积分时间相当于carrierloop_fll_cohtime/2
            if(mscountL1_bitperiod<=bitperiod/2)
                corrcohsum_carr_FLL_first=corrcohsum_carr_FLL_first+corrsum(mid+0); %为FLL作更新  每一个电文周期的前一半
            else
                corrcohsum_carr_FLL_second=corrcohsum_carr_FLL_second+corrsum(mid+0);%为FLL作更新  每一个电文周期的后一半
            end
            corrsum = zeros(1,LengthCorreArray);                   %每ms累加值清零
            if(mode == 0)
                %PLL
                if( mod(mscountL1_bitperiod,carrierloop_cohtime) == 0)  %载波环路更新
                    P = corrcohsum_carr(mid);                     %使用P路
                    pd = atan(imag(P)/real(P))/2/pi;              %鉴相结果
                    indexpd=floor(loopcnt/carrierloop_cohtime);
                    PD(indexpd)= pd;
                    delta_pd_n1 = pd - pd_1;                     %二阶锁相环使用
                    delta_pd_n2 = pd_1 - pd_2;                   %三阶锁相环使用
                    delta_pd_n3 = pd_2 - pd_3;
                    if(flag == 1)
                        if(abs(delta_pd_n1) >= 0.25)
                            if(delta_pd_n1 > 0)
                                delta_pd_n1 = (delta_pd_n1-0.5);
                            else
                                delta_pd_n1 = (delta_pd_n1+0.5);
                            end
                        end
                        if(abs(delta_pd_n2) >= 0.25)
                            if(delta_pd_n2 > 0)
                                delta_pd_n2 = (delta_pd_n2-0.5);
                            else
                                delta_pd_n2 = (delta_pd_n2+0.5);
                            end
                        end
                        
                        if(abs(delta_pd_n3) >= 0.25)
                            if(delta_pd_n3 > 0)
                                delta_pd_n3 = (delta_pd_n3-0.5);
                            else
                                delta_pd_n3 = (delta_pd_n3+0.5);
                            end
                        end
                        
                    end
                    if(carrierlooporder == 2)
                        vn = vn_1 + cp1_2nd * (delta_pd_n1) + cp2_2nd*T*carrierloop_cohtime/2*(pd + pd_1);
                    elseif(carrierlooporder == 3)
                        vn = 2*vn_1 - vn_2 + cp1_3rd * (delta_pd_n1-delta_pd_n2) + ...
                            cp2_3rd *T*carrierloop_cohtime/2* (delta_pd_n1+delta_pd_n2)	+ cp3_3rd *(T*carrierloop_cohtime/2)^2* (pd + 2*pd_1 + pd_2);
                    elseif(carrierlooporder == 4)
                        vn = 3*vn_1 - 3*vn_2 + vn_3  + cp1_4th * (delta_pd_n1-2*delta_pd_n2+delta_pd_n3 ) + ...%pd-3pd_1+3pd_2-pd_3
                            cp2_4th *(T*carrierloop_cohtime/2)* (delta_pd_n1-delta_pd_n3)	+ ...                                           %pd-pd_1-pd_2+pd_3
                            cp3_4th *(T*carrierloop_cohtime/2)^2* (delta_pd_n1+delta_pd_n3+2*delta_pd_n2 )+ ...                             %pd+pd_1-pd_2-pd_3
                            cp4_4th *(T*carrierloop_cohtime/2)^3* (pd+3*pd_1+3*pd_2+pd_3);
                        
                    end
                    pll_lp_out = vn - vn_1;
                    carrier_dco_deltaL = carrier_dco_deltaL + pll_lp_out/Freq_sample;
                    
                    %             carrier_dco_deltaL*Freq_sample-Fc
                    vn_3 = vn_2;
                    vn_2 = vn_1;
                    vn_1 = vn;
                    pd_3 = pd_2;
                    pd_2 = pd_1;
                    pd_1 = pd;
                    corrcohsum_carr = zeros(1,LengthCorreArray);
                    
                end
                %FLL
            elseif(mode==1)
                if( mod(mscountL1_bitperiod,carrierloop_fll_cohtime) == 0)
                    %FLL鉴相 使用二象限反正切 每个比特周期内算一次dot与CRoss
                    carrierloop_fll_noncohcnt = carrierloop_fll_noncohcnt + 1;
                    cross_dot = cross_dot + conj(corrcohsum_carr_FLL_first)*corrcohsum_carr_FLL_second;
                    if(carrierloop_fll_noncohcnt >= carrierloop_fll_noncohtime)
                        
                        
                        fd = atan2(imag(cross_dot),real(cross_dot))/(T*carrierloop_fll_real_cohtime*2*pi);%使用4象限反正切
                        %                     fd=carrierloop_fll_noncohtime*fd;
                        cross_dot = 0;
                        carrierloop_fll_noncohcnt = 0;
                        indexPDFLL=floor(loopcnt/(carrierloop_fll_noncohtime*carrierloop_fll_cohtime));
                        PDFLL(indexPDFLL)= fd;
                        %根据阶数改变NCO
                        if(carrierloopFLLorder == 1)
                            fo_n = fo_n_1 + cf1_1st * (fd + fd_1);
                        elseif(carrierloopFLLorder == 2)
                            fo_n=   2*fo_n_1-fo_n_2+cf1_2nd*(fd-fd_2)+cf2_2nd*(fd+2*fd_1+fd_2);
                        elseif(carrierloopFLLorder == 3)
                            fo_n = 3*fo_n_1 - 3*fo_n_2 + fo_n_3+ cf1_3rd * (fd - fd_1 - fd_2 + fd_3)+ cf2_3rd * (fd + fd_1 - fd_2 - fd_3)+ cf3_3rd * (fd + 3*fd_1 + 3*fd_2 + fd_3);
                        end
                        carrier_dco_deltaL = carrier_dco_deltaL + (fo_n-fo_n_1)/Freq_sample;
                        
                        fd_3 = fd_2;
                        fd_2 = fd_1;
                        fd_1 = fd;
                        fo_n_3 = fo_n_2;
                        fo_n_2 = fo_n_1;
                        fo_n_1 = fo_n;
                        
                    end
                    corrcohsum_carr_FLL_first = 0;
                    corrcohsum_carr_FLL_second = 0;
                end
            elseif(mode==2)
                if(mod(mstime,2)==1)
                    %% ———— 新增：mode==2 卡尔曼滤波分支 ————
                    
                    % —— 第一次进入时做一次初始化 ——
                    if ~exist('kf_initialized','var') || ~kf_initialized
                        % 时间步长（同环路调用率）
                        T_kf = 2*1e-3;
                        T_coh=1*1e-3;
                        %xk=Φ*xk-1+B*G
                        % 状态转移矩阵 Φ (4×4)
                        Phi_kf = [1, T_kf, T_kf^2/2, T_kf^3/6;
                            0,     1,      T_kf,    T_kf^2/2;
                            0,     0,        1,         T_kf;
                            0,     0,        0,           1];
                        %
                        B_kf=[-1 ,-T_kf;
                            0,0 ;
                            0,0;
                            0,0;];
                        %G外加控制量 有deltaagnel相位差调整量==Xk-1(1),以及f_nco（k-1）
                        f_nco=LocalDop;%f_nco单位是频率
                        deltaagnel=0;
                        % 过程噪声强度 σ?： A_max是jerk PSD 单位 周期^2/s^5
                        A_max = 50*100;
                        sigma_a = (A_max) * T_kf ;
                        Q_kf = sigma_a * [T_kf^6/252, T_kf^5/72, T_kf^4/30, T_kf^3/24;
                            T_kf^5/72, T_kf^4/20,  T_kf^3/8, T_kf^2/6;
                            T_kf^4/30,  T_kf^3/8,    T_kf^2/3,    T_kf/2;
                            T_kf^3/24,  T_kf^2/6,    T_kf/2,       1];
                        
                        % 观测矩阵 H
                        H_kf = [1, T_kf*1/2, T_kf^2/6, T_kf^3/24];
                        
                        % 测量噪声 R
                        COHNUM=Freq_sample*T_coh;
                        R_kf=1/(COHNUM)  ;
                        x_kf = zeros(4,1);%
                        x_kf(2) = LocalDop;     % 当前多普勒频率估计
                        P_kf = diag([5, 20^2, (20/2)^2, 10]);% P_kf = diag([pi^2, 100^2, 50^2, 10^2]);%
                        kf_initialized = true;
                        
                    end
                    
                    % —— 采样到新的相位鉴相器输出 pd ——
                    P = corrcohsum_carr(mid);
                    %                 P=P/AMP;
                    pd = atan(imag(P)/real(P)) / (2*pi);  %单位输出是 周期数
                    %计算索引
                    indexLKFpd=floor(loopcnt/2);
                    
                    % —— 1. 预测 ——
                    x_pred = Phi_kf * x_kf+B_kf*[deltaagnel;f_nco];
                    P_pred = Phi_kf * P_kf * Phi_kf' + Q_kf;
                    %计算观测的残差 inv部分，也就是zk-hk*xpred-ck*uk  单位是周期数
                    rawinv=H_kf*x_pred-T_kf/2*(f_nco);
                    innov=pd-rawinv;
                    PD(indexLKFpd)= pd;
                    % —— 2. 更新 ——
                    K_kf = P_pred * H_kf' / (H_kf * P_pred * H_kf' + R_kf);
                    x_kf = x_pred + K_kf * innov;
                    P_kf = (eye(4) - K_kf * H_kf) * P_pred;
                    f_nco = x_kf(2) + x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    carrier_dco_deltaL = (Fc + f_nco) / Freq_sample;
                    deltaagnel= x_kf(1);
                    carrier_dco_phaseL = carrier_dco_phaseL+deltaagnel;%相位调整量，直接加上本次滤波之后的相位差量，必须！
                    % 保持其余状态变量推进
                    corrcohsum_carr = zeros(1,LengthCorreArray);
                    
                    fid = fopen(filename,'a+');
                    fprintf(fid,'loopcnt %06d carr %+010.2f  pd%.3f  innov%.3f\n',...
                        loopcnt,carrier_dco_deltaL*Freq_sample-Fc,pd,innov);
                    fclose(fid);
                    %                 fprintf('P_k对角线: [%.2e, %.2e, %.2e, %.2e]\n', diag(P_kf));
                    
                    
                end
            elseif(mode==3)
                if(mod(mstime,2)==1)
                    %% ———— 新增：mode==3 EKF滤波分支 ————
                    
                    % —— 第一次进入时做一次初始化 ——
                    if ~exist('Ekf_initialized','var') || ~Ekf_initialized
                        % 时间步长（同环路调用率）
                        T_kf = 2*1e-3;
                        %xk=Φ*xk-1+B*G
                        % 状态转移矩阵 Φ (4×4)
                        Phi_kf = [1, T_kf, T_kf^2/2, T_kf^3/6;
                            0,     1,      T_kf,    T_kf^2/2;
                            0,     0,        1,         T_kf;
                            0,     0,        0,           1];
                        %
                        B_kf=[-1 ,-T_kf;
                            0,0 ;
                            0,0;
                            0,0;];
                        %G外加控制量 有deltaagnel相位差调整量==Xk-1(1),以及f_nco（k-1）
                        f_nco=LocalDop;%f_nco单位是频率
                        deltaagnel=0;
                        % 过程噪声强度 σ?：可以根据最大加加速度 A_max 来估计
                        A_max = 5000;
                        sigma_a = (A_max) * T_kf / 2;
                        Q_kf = sigma_a * [T_kf^6/252, T_kf^5/72, T_kf^4/30, T_kf^3/24;
                            T_kf^5/72, T_kf^4/20,  T_kf^3/8, T_kf^2/6;
                            T_kf^4/30,  T_kf^3/8,    T_kf^2/3,    T_kf/2;
                            T_kf^3/24,  T_kf^2/6,    T_kf/2,       1];
                        
                        % 观测矩阵 H  与KF不同，观测矩阵每一次都是会改变的观测方程对各个状态变量的偏导数
                        H_kf = zeros(2,4);
                        %                     HA=A(AcqSatNUM)*T_kf*Freq_sample;% AT
                        %H11=AT*sin(2piθe k|k-1 )*(-2pi)   H21=AT*cos(2piθe k|k-1 )*(+2pi)
                        % 测量噪声 R
                        %R_kf=simga0^2/Ncoh=simga0^2/(Tcoh*fs)=B/(Tcoh*fs*fs)
                        %                     rtemp=4.092e6/(2*1e-3*Freq_sample^2);
                        rtemp=1/(Freq_sample*T_kf);
                        R_kf=[ rtemp,0;0, rtemp];
                        % 状态初始化 x = [θe; f; f?; f?]，从当前环路值出发
                        x_kf = zeros(4,1);
                        %                     x_kf(1) = mod(LocalDop*T_kf,1);            % 残余相位初值
                        x_kf(2) = LocalDop;     % 当前多普勒频率估计
                        % f?, f? 初始设 0
                        P_kf = diag([5, 20^2, (20/2)^2, 10]);% P_kf = diag([pi^2, 100^2, 50^2, 10^2]);
                        
                        Ekf_initialized = true;
                        A1=A(AcqSatNUM);
                    end
                    
                    % —— 采样到新的相位鉴相器输出 pd ——
                    P = corrcohsum_carr(mid);
                    %电文的影响怎么办 不能简单的取绝对值    innov=Z_kf-H_forinnov;  一定会影响
                    %剥离？？？
                    I_k=real(P)/(Freq_sample*T_kf/2);
                    Q_k=imag(P)/(Freq_sample*T_kf/2);
                    
                    
                    % —— 1. 预测 ——
                    x_pred = Phi_kf * x_kf+B_kf*[deltaagnel;f_nco];
                    P_pred = Phi_kf * P_kf * Phi_kf' + Q_kf;
                    %计算观测的残差 inv部分，也就是zk-hk*xpred-ck*uk  单位是周期数
                    Z_kf=[I_k;Q_k];
                    %H11=AT*sin(2piθe k|k-1 )*(-2pi)
                    H_kf(1,1)=A1*sin(2*pi*x_pred(1))*(-2*pi);
                    H_kf(2,1)=A1*cos(2*pi*x_pred(1))*(2*pi);
                    
                    H_kf(1,2)=A1*sin(2*pi*x_pred(2))*(-pi)*T_kf;
                    H_kf(2,2)=A1*cos(2*pi*x_pred(2))*(pi)*T_kf;
                    
                    H_kf(1,3)=A1*sin(2*pi*x_pred(3))*(-pi)*T_kf*T_kf/3;
                    H_kf(2,3)=A1*cos(2*pi*x_pred(3))*(pi)*T_kf*T_kf/3;
                    
                    H_kf(1,4)=A1*sin(2*pi*x_pred(4))*(-pi)*T_kf*T_kf*T_kf/12;
                    H_kf(2,4)=A1*cos(2*pi*x_pred(4))*(pi)*T_kf*T_kf*T_kf/12;
                    %校正
                    H_forinnov=[  A1*cos(2 * pi * x_pred(1)); A1*sin(2 * pi * x_pred(1))];
                    innov=Z_kf-H_forinnov;
                    
                    % —— 2. 更新 ——
                    K_kf = P_pred * H_kf' / (H_kf * P_pred * H_kf' + R_kf);
                    x_kf = x_pred + K_kf * innov;
                    P_kf = (eye(4) - K_kf * H_kf) * P_pred;
                    f_nco = x_kf(2) + x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    carrier_dco_deltaL = (Fc + f_nco) / Freq_sample;
                    deltaagnel= x_kf(1);
                    carrier_dco_phaseL = carrier_dco_phaseL+deltaagnel;%相位调整量，直接加上本次滤波之后的相位差量，必须！
                    % 保持其余状态变量推进
                    corrcohsum_carr = zeros(1,LengthCorreArray);
                    
                    
                    fid = fopen(filename,'a+');
                    fprintf(fid,'loopcnt %06d carr %+010.2f  pd%.3f\n',...
                        loopcnt,carrier_dco_deltaL*Freq_sample-Fc,pd);
                    fclose(fid);
                    %                 fprintf('P_k对角线: [%.2e, %.2e, %.2e, %.2e]\n', diag(P_kf));
                    
                end
            elseif(mode==4)%不剥离电文  暂时不行
                if(mod(mstime,2)==1)
                    %% ———— 新增：mode==4 EKF滤波分支不剥离电文增加电文的预测 ————
                    
                    % —— 第一次进入时做一次初始化 ——
                    if ~exist('Ekf_initializedbit','var') || ~Ekf_initializedbit
                        % 时间步长（同环路调用率）
                        T_kf = 2*1e-3;
                        %xk=Φ*xk-1+B*G
                        % 状态转移矩阵 Φ (4×4)
                        Phi_kf = [1, T_kf, T_kf^2/2, T_kf^3/6,0;
                            0,     1,      T_kf,    T_kf^2/2,0;
                            0,     0,        1,         T_kf,0;
                            0,     0,        0,           1,0;
                            0,0,0,0,0.95];
                        %
                        B_kf=[-1 ,-T_kf;
                            0,0 ;
                            0,0;
                            0,0;
                            0,0;];
                        %G外加控制量 有deltaagnel相位差调整量==Xk-1(1),以及f_nco（k-1）
                        f_nco=LocalDop;%f_nco单位是频率
                        deltaagnel=0;
                        % 过程噪声强度 σ?：可以根据最大加加速度 A_max 来估计
                        A_max = 1000;
                        sigma_a = (A_max) * T_kf / 2;
                        Q_kf = sigma_a * [T_kf^6/252, T_kf^5/72, T_kf^4/30, T_kf^3/24,0;
                            T_kf^5/72, T_kf^4/20,  T_kf^3/8, T_kf^2/6,0;
                            T_kf^4/30,  T_kf^3/8,    T_kf^2/3,    T_kf/2,0;
                            T_kf^3/24,  T_kf^2/6,    T_kf/2,       1,0;
                            0 ,0 ,0 ,0 ,1;                               ];
                        
                        % 观测矩阵 H  与KF不同，观测矩阵每一次都是会改变的观测方程对各个状态变量的偏导数
                        H_kf = zeros(2,5);
                        HA=A(AcqSatNUM)*T_kf*Freq_sample;% AT
                        %H11=AT*sin(2piθe k|k-1 )*(-2pi)   H21=AT*cos(2piθe k|k-1 )*(+2pi)
                        % 测量噪声 R
                        %R_kf=simga0^2/Ncoh=simga0^2/(Tcoh*fs)=B/(Tcoh*fs*fs)
                        rtemp=4.092e6/(2*1e-3*Freq_sample^2);
                        R_kf=[ rtemp,0;0, rtemp];
                        % 状态初始化 x = [θe; f; f?; f?,Dk]，从当前环路值出发
                        x_kf = zeros(5,1);
                        %                     x_kf(1) = mod(LocalDop*T_kf,1);            % 残余相位初值
                        x_kf(2) = LocalDop;     % 当前多普勒频率估计
                        % f?, f? 初始设 0
                        x_kf(5)=1;
                        P_kf = diag([5, 20^2, (20/2)^2, 10,2]);% P_kf = diag([pi^2, 100^2, 50^2, 10^2]);
                        predictbit=[];%存放预测的比特看看对不对？
                        Ekf_initializedbit = true;
                        
                    end
                    
                    % —— 采样到新的相位鉴相器输出 pd ——
                    P = corrcohsum_carr(mid);
                    %电文的影响怎么办 不能简单的取绝对值    innov=Z_kf-H_forinnov;  一定会影响
                    %剥离？？？
                    I_k=2*real(P)/HA;
                    Q_k=2*imag(P)/HA;
                    
                    
                    % —— 1. 预测 ——
                    x_pred = Phi_kf * x_kf+B_kf*[deltaagnel;f_nco];
                    P_pred = Phi_kf * P_kf * Phi_kf' + Q_kf;
                    %计算观测的残差 inv部分，也就是zk-hk*xpred-ck*uk  单位是周期数
                    Z_kf=[I_k;Q_k];
                    %H11=AT*sin(2piθe k|k-1 )*(-2pi)
                    H_kf(1,1)=sin(2*pi*x_pred(1))*(-2*pi)*x_pred(5);
                    H_kf(2,1)=cos(2*pi*x_pred(1))*(2*pi)*x_pred(5);
                    
                    %                 H_kf(1,2)=sin(2*pi*x_pred(2))*(-pi)*T_kf;
                    %                 H_kf(2,2)=cos(2*pi*x_pred(2))*(pi)*T_kf;
                    %
                    %                 H_kf(1,3)=sin(2*pi*x_pred(3))*(-pi)*T_kf*T_kf/3;
                    %                 H_kf(2,3)=cos(2*pi*x_pred(3))*(pi)*T_kf*T_kf/3;
                    %
                    %                 H_kf(1,4)=sin(2*pi*x_pred(4))*(-pi)*T_kf*T_kf*T_kf/12;
                    %                 H_kf(2,4)=cos(2*pi*x_pred(4))*(pi)*T_kf*T_kf*T_kf/12;
                    
                    H_kf(1,5)=cos(2*pi*x_pred(5));
                    H_kf(2,5)=sin(2*pi*x_pred(5));
                    
                    %校正
                    H_forinnov=[  x_pred(5)*cos(2 * pi * x_pred(1)); x_pred(5)*sin(2 * pi * x_pred(1))];
                    innov=Z_kf-H_forinnov;
                    
                    % —— 2. 更新 ——
                    K_kf = P_pred * H_kf' / (H_kf * P_pred * H_kf' + R_kf);
                    x_kf = x_pred + K_kf * innov;
                    x_kf(5)=sign(x_kf(5));
                    %                 predictbit=[predictbit x_kf(5)];
                    P_kf = (eye(5) - K_kf * H_kf) * P_pred;
                    f_nco = x_kf(2) + x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    carrier_dco_deltaL = (Fc + f_nco) / Freq_sample;
                    deltaagnel= x_kf(1);
                    carrier_dco_phaseL = carrier_dco_phaseL+deltaagnel;%相位调整量，直接加上本次滤波之后的相位差量，必须！
                    % 保持其余状态变量推进
                    corrcohsum_carr = zeros(1,LengthCorreArray);
                    
                    
                    fid = fopen(filename,'a+');
                    fprintf(fid,'loopcnt %06d carr %+010.2f  pd%.3f\n',...
                        loopcnt,carrier_dco_deltaL*Freq_sample-Fc,pd);
                    fclose(fid);
                    %                 fprintf('P_k对角线: [%.2e, %.2e, %.2e, %.2e]\n', diag(P_kf));
                    
                end
                
            elseif(mode==5)%CSK解调，不调整
                %                 FCARR = [FCARR carrier_dco_deltaL*Freq_sample-Fc];
                %                 FCODE = [FCODE;loopcnt code_dco_deltaL*Freq_sample-Freq_code];
            elseif(mode==6)
                if(mod(mstime,2)==1)
                    %% ———— 新增：mode==6 UKF滤波分支剥离电文 ————
                    
                    % —— 第一次进入时做一次初始化 ——
                    if ~exist('Ukf_initializedbit','var') || ~Ukf_initializedbit
                        % 时间步长（同环路调用率）
                        T_kf = 2*1e-3;
                        % —— UKF 参数 ——
                        nukf = 4;                        % 状态维数
                        alpha = 1e-3;  kappa = 0;  beta = 2;
                        lambda = alpha^2*(nukf+kappa)-nukf;
                        wm = [lambda/(nukf+lambda) repmat(1/(2*(nukf+lambda)),1,2*nukf)];
                        wc = wm;
                        wc(1) = wm(1)+(1-alpha^2+beta);
                        
                        Phi_kf = [1, T_kf, T_kf^2/2, T_kf^3/6;
                            0,     1,      T_kf,    T_kf^2/2;
                            0,     0,        1,         T_kf;
                            0,     0,        0,           1];
                        %
                        B_kf=[-1 ,-T_kf;
                            0,0 ;
                            0,0;
                            0,0;];
                        %G外加控制量 有deltaagnel相位差调整量==Xk-1(1),以及f_nco（k-1）
                        f_nco=LocalDop;%f_nco单位是频率
                        
                        % 过程噪声强度 σ?：可以根据最大加加速度 A_max 来估计
                        A_max = 500;
                        sigma_a = (A_max) * T_kf / 2;
                        Q_kf = sigma_a * [T_kf^6/252, T_kf^5/72, T_kf^4/30, T_kf^3/24;
                            T_kf^5/72, T_kf^4/20,  T_kf^3/8, T_kf^2/6;
                            T_kf^4/30,  T_kf^3/8,    T_kf^2/3,    T_kf/2;
                            T_kf^3/24,  T_kf^2/6,    T_kf/2,       1;  ];
                        
                        % 观测矩阵 H  与KF不同，观测矩阵每一次都是会改变的观测方程对各个状态变量的偏导数
                        H_kf = zeros(2,4);
                        HA=A(AcqSatNUM)*T_kf*Freq_sample;% AT
                        deltaagnel=0;
                        % 测量噪声 R
                        %R_kf=simga0^2/Ncoh=simga0^2/(Tcoh*fs)=B/(Tcoh*fs*fs)
                        rtemp=4.092e6/(2*1e-3*Freq_sample^2);
                        R_kf=[ rtemp,0;0, rtemp];
                        % 状态初始化 x = [θe; f; f?; f?,Dk]，从当前环路值出发
                        x_kf = zeros(4,1);
                        %                     x_kf(1) = mod(LocalDop*T_kf,1);            % 残余相位初值
                        x_kf(2) = LocalDop;     % 当前多普勒频率估计
                        % f?, f? 初始设 0
                        
                        P_kf = diag([5, 20^2, (20/2)^2, 10]);% P_kf = diag([pi^2, 100^2, 50^2, 10^2]);
                        
                        Ukf_initializedbit = true;
                        % —— 初始状态和协方差 ——
                        x_ukf = x_kf;                 % 从之前分支继承或自定义
                        P_ukf = P_kf;
                        Q_ukf = Q_kf;  R_ukf = R_kf;  % I/Q 噪声协方差
                    end
                    % —— 生成 2n+1 sigma 点矩阵 —— 相当于 X （K-1） 从k-1最优估计选取 2n+1 个生成sigma
                    S = chol((nukf+lambda)*P_ukf,'lower');
                    sigmaukf = [ x_ukf, x_ukf+S, x_ukf-S ];  % 大小 n×(2n+1)
                    
                    % —— UKF 预测 ——  相当于 XK|k-1
                    % 利用状态转移矩阵将2n+1个生成sigma进行预测，预测到第k个时刻的 2n+1 个预测sigma点
                    for i=1:2*nukf+1
                        xi = sigmaukf(:,i);
                        % 状态转移：同 LKF 的 Phi*xi + B*u
                        sigma_pred(:,i) = Phi_kf*xi + B_kf*[deltaagnel;f_nco];
                    end
                    %预测sigma点计算均值与协方差
                    x_pred = sigma_pred * wm';              % 均值预测
                    P_pred = Q_ukf;
                    for i=1:2*nukf+1
                        d = sigma_pred(:,i)-x_pred;
                        P_pred = P_pred + wc(i)*(d*d');
                    end
                    
                    % —— UKF 观测预测 ——
                    %预测 2n+1 sigma点带入观测方程中，得到预测的观测值
                    for i=1:2*nukf+1
                        th = sigma_pred(1,i);
                        f  = sigma_pred(2,i);
                        % 非线性观测 h(x) h与EKF是一样的，都是I Q观测值，注意电文
                        z_sigma(:,i) = [cos(2*pi*th); sin(2*pi*th)];
                    end
                    z_pred = z_sigma * wm';                 % 观测均值
                    S_ukf = R_ukf;
                    Pxz = zeros(nukf,2);
                    for i=1:2*nukf+1
                        dz = z_sigma(:,i)-z_pred;
                        dx = sigma_pred(:,i)-x_pred;
                        S_ukf  = S_ukf + wc(i)*(dz*dz');
                        Pxz    = Pxz   + wc(i)*(dx*dz');
                    end
                    
                    % —— UKF 更新 ——
                    K_ukf = Pxz / S_ukf;
                    P_ukf = P_pred - K_ukf*S_ukf*K_ukf';
                    % 获取真实测量 I/Q
                    %注意电文！！
                    Pxy = corrcohsum_carr(mid);
                    z    = [2*real(Pxy)/HA; 2*imag(Pxy)/HA];
                    innov = z - z_pred;
                    x_ukf = x_pred + K_ukf * innov;
                    
                    
                    f_nco = x_ukf(2) + x_ukf(3)*T_kf + x_ukf(4)*T_kf^2/2;
                    carrier_dco_deltaL = (Fc + f_nco) / Freq_sample;
                    deltaagnel= x_ukf(1);
                    carrier_dco_phaseL = carrier_dco_phaseL+deltaagnel;%相位调整量，直接加上本次滤波之后的相位差量，必须！
                    % 保持其余状态变量推进
                    corrcohsum_carr = zeros(1,LengthCorreArray);
                    
                    
                    fid = fopen(filename,'a+');
                    fprintf(fid,'loopcnt %06d carr %+010.2f  pd%.3f\n',...
                        loopcnt,carrier_dco_deltaL*Freq_sample-Fc,pd);
                    fclose(fid);
                    %                 fprintf('P_k对角线: [%.2e, %.2e, %.2e, %.2e]\n', diag(P_kf));
                    
                end
                
                
            elseif(mode==7)%预留
                
            elseif(mode==8)%LKF 1ms 预测，2ms矫正
                
                %% ———— 新增：mode==8 卡尔曼滤波分支 ———— 1ms预测，2ms矫正
                
                % —— 第一次进入时做一次初始化 ——
                if ~exist('kf_initialized','var') || ~kf_initialized
                    % 时间步长（同环路调用率）
                    T_kf = 1*1e-3;
                    T_coh=T_kf;
                    %xk=Φ*xk-1+B*G
                    % 状态转移矩阵 Φ (4×4)
                    Phi_kf = [1, T_kf, T_kf^2/2, T_kf^3/6;
                        0,     1,      T_kf,    T_kf^2/2;
                        0,     0,        1,         T_kf;
                        0,     0,        0,           1];
                    %
                    B_kf=[-1 ,-T_kf;
                        0,0 ;
                        0,0;
                        0,0;];
                    %G外加控制量 有deltaagnel相位差调整量==Xk-1(1),以及f_nco（k-1）
                    f_nco=LocalDop;%f_nco单位是频率
                    deltaagnel=0;
                    % 过程噪声强度 σ?：可以根据最大加加速度 A_max 来估计
                    A_max = 500;
                    sigma_a = (A_max) * T_kf / 2;
                    Q_kf = sigma_a * [T_kf^6/252, T_kf^5/72, T_kf^4/30, T_kf^3/24;
                        T_kf^5/72, T_kf^4/20,  T_kf^3/8, T_kf^2/6;
                        T_kf^4/30,  T_kf^3/8,    T_kf^2/3,    T_kf/2;
                        T_kf^3/24,  T_kf^2/6,    T_kf/2,       1];
                    
                    % 观测矩阵 H
                    T_kf_for_obv=2e-3/2;
                    H_kf = [1, T_kf_for_obv*1/2, T_kf_for_obv^2/6, T_kf_for_obv^3/24];
                    
                    %                 % 测量噪声 R
                    COHNUM=Freq_sample*T_coh;
                    R_kf=1/COHNUM;
                    % 状态初始化 x = [θe; f; f?; f?]，从当前环路值出发
                    x_kf = zeros(4,1);
                    %                     x_kf(1) = mod(LocalDop*T_kf,1);            % 残余相位初值
                    x_kf(2) = LocalDop;     % 当前多普勒频率估计
                    % f?, f? 初始设 0
                    P_kf = diag([5, 20^2, (20/2)^2, 10]);% P_kf = diag([pi^2, 100^2, 50^2, 10^2]);
                    f_csk=LocalDop;
                    kf_initialized = true;
                    
                end
                
                % —— 采样到新的相位鉴相器输出 pd ——
                P = corrcohsum_carr(mid);
                pd = atan(imag(P)/real(P)) / (2*pi);  %单位输出是 周期数
                
                
                
                % —— 1. 预测 ——
                x_pred = Phi_kf * x_kf+B_kf*[deltaagnel;f_nco];
                P_pred = Phi_kf * P_kf * Phi_kf' + Q_kf;
                %计算观测的残差 inv部分，也就是zk-hk*xpred-ck*uk  单位是周期数
                if(mod(mstime,2)==1)
                    rawinv=H_kf*x_pred-T_kf/2*(f_nco);
                    innov=pd-rawinv;
                    K_kf = P_pred * H_kf' / (H_kf * P_pred * H_kf' + R_kf);
                    x_kf = x_pred + K_kf * innov;
                    P_kf = (eye(4) - K_kf * H_kf) * P_pred;
                    temp=[0 0 T_kf 1/2*T_kf*T_kf];
                    tuo=temp*P_kf*temp'+566*0.001*5;
                    
                    f_nco = x_kf(2) + x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    
                    if((f_nco-f_csk)>50)
                        f_nco=f_csk+50;
                    elseif((f_nco-f_csk)<-50)
                        f_nco=f_csk-50;
                    end
                    f_nco_bpsk=f_nco;
                else %外推算法限制条件
                    fpred= x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    flim=max(min(fpred,tuo),-tuo);
                    f_nco=f_nco_bpsk+flim;
                    f_csk=f_nco;
                end
                
                % —— 2. 更新 ——
                
                
                carrier_dco_deltaL = (Fc + f_nco) / Freq_sample;
                deltaagnel= x_kf(1);
                carrier_dco_phaseL = carrier_dco_phaseL+deltaagnel;%相位调整量，直接加上本次滤波之后的相位差量，必须！
                % 保持其余状态变量推进
                corrcohsum_carr = zeros(1,LengthCorreArray);
                
                
                fid = fopen(filename,'a+');
                %             fprintf(fid,'loopcnt %06d carr %+010.2f  pd%.3f  innov%.3f\n',...
                %                 loopcnt,carrier_dco_deltaL*Freq_sample-Fc,pd,innov);
                fclose(fid);
                %                 fprintf('P_k对角线: [%.2e, %.2e, %.2e, %.2e]\n', diag(P_kf));
                
                
                
            elseif(mode==9)
                
                %% ———— 新增：mode==9 EKF滤波分支 ————
                %1ms预测 2ms校正
                % —— 第一次进入时做一次初始化 ——
                if ~exist('Ekf_initialized','var') || ~Ekf_initialized
                    % 时间步长（同环路调用率）
                    T_kf = 1*1e-3;
                    %xk=Φ*xk-1+B*G
                    % 状态转移矩阵 Φ (4×4)
                    Phi_kf = [1, T_kf, T_kf^2/2, T_kf^3/6;
                        0,     1,      T_kf,    T_kf^2/2;
                        0,     0,        1,         T_kf;
                        0,     0,        0,           1];
                    %
                    B_kf=[-1 ,-T_kf;
                        0,0 ;
                        0,0;
                        0,0;];
                    %G外加控制量 有deltaagnel相位差调整量==Xk-1(1),以及f_nco（k-1）
                    f_nco=LocalDop;%f_nco单位是频率
                    deltaagnel=0;
                    % 过程噪声强度 σ?：可以根据最大加加速度 A_max 来估计
                    A_max = 5000;
                    sigma_a = (A_max) * T_kf / 2;
                    Q_kf = sigma_a * [T_kf^6/252, T_kf^5/72, T_kf^4/30, T_kf^3/24;
                        T_kf^5/72, T_kf^4/20,  T_kf^3/8, T_kf^2/6;
                        T_kf^4/30,  T_kf^3/8,    T_kf^2/3,    T_kf/2;
                        T_kf^3/24,  T_kf^2/6,    T_kf/2,       1];
                    
                    % 观测矩阵 H  与KF不同，观测矩阵每一次都是会改变的观测方程对各个状态变量的偏导数
                    H_kf = zeros(2,4);
                    rtemp=1/(Freq_sample*T_kf);
                    R_kf=[ rtemp,0;0, rtemp];
                    % 状态初始化 x = [θe; f; f?; f?]，从当前环路值出发
                    x_kf = zeros(4,1);
                    %                     x_kf(1) = mod(LocalDop*T_kf,1);            % 残余相位初值
                    x_kf(2) = LocalDop;     % 当前多普勒频率估计
                    
                    P_kf = diag([5, 20^2, (20/2)^2, 10]);% P_kf = diag([pi^2, 100^2, 50^2, 10^2]);
                    f_csk=LocalDop;
                    Ekf_initialized = true;
                    A1=A(AcqSatNUM);
                    
                end
                
                % —— 采样到新的相位鉴相器输出 pd ——
                P = corrcohsum_carr(mid);
                %电文的影响怎么办 不能简单的取绝对值    innov=Z_kf-H_forinnov;  一定会影响
                %剥离
                I_k=real(P)/(Freq_sample*T_kf);
                Q_k=imag(P)/(Freq_sample*T_kf);
                
                % —— 1. 预测 ——
                x_pred = Phi_kf * x_kf+B_kf*[deltaagnel;f_nco];
                P_pred = Phi_kf * P_kf * Phi_kf' + Q_kf;
                %计算观测的残差 inv部分，也就是zk-hk*xpred-ck*uk  单位是周期数
                if(mod(mstime,2)==1)
                    Z_kf=[I_k;Q_k];
                    %H11=AT*sin(2piθe k|k-1 )*(-2pi)
                    H_kf(1,1)=A1*sin(2*pi*x_pred(1))*(-2*pi);
                    H_kf(2,1)=A1*cos(2*pi*x_pred(1))*(2*pi);
                    
                    H_kf(1,2)=A1*sin(2*pi*x_pred(2))*(-pi)*T_kf;
                    H_kf(2,2)=A1*cos(2*pi*x_pred(2))*(pi)*T_kf;
                    
                    H_kf(1,3)=A1*sin(2*pi*x_pred(3))*(-pi)*T_kf*T_kf/3;
                    H_kf(2,3)=A1*cos(2*pi*x_pred(3))*(pi)*T_kf*T_kf/3;
                    
                    H_kf(1,4)=A1*sin(2*pi*x_pred(4))*(-pi)*T_kf*T_kf*T_kf/12;
                    H_kf(2,4)=A1*cos(2*pi*x_pred(4))*(pi)*T_kf*T_kf*T_kf/12;
                    %校正
                    H_forinnov=[  A1*cos(2 * pi * x_pred(1)); A1*sin(2 * pi * x_pred(1))];
                    innov=Z_kf-H_forinnov;
                    
                    % —— 2. 更新 ——
                    K_kf = P_pred * H_kf' / (H_kf * P_pred * H_kf' + R_kf);
                    x_kf = x_pred + K_kf * innov;
                    P_kf = (eye(4) - K_kf * H_kf) * P_pred;
                    
                    temp=[0 0 T_kf 1/2*T_kf*T_kf];
                    tuo=temp*P_kf*temp'+566*0.001*5;
                    f_nco = x_kf(2) + x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    
                    if((f_nco-f_csk)>50)
                        f_nco=f_csk+50;
                    elseif((f_nco-f_csk)<-50)
                        f_nco=f_csk-50;
                    end
                    f_nco_bpsk=f_nco;
                    
                else
                    fpred= x_kf(3)*T_kf + x_kf(4)*T_kf^2/2;
                    flim=max(min(fpred,tuo),-tuo);
                    f_nco=f_nco_bpsk+flim;
                    f_csk=f_nco;
                end
                
                carrier_dco_deltaL = (Fc + f_nco) / Freq_sample;
                deltaagnel= x_kf(1);
                carrier_dco_phaseL = carrier_dco_phaseL+deltaagnel;%相位调整量，直接加上本次滤波之后的相位差量，必须！
                % 保持其余状态变量推进
                corrcohsum_carr = zeros(1,LengthCorreArray);
                
                
                fid = fopen(filename,'a+');
                fprintf(fid,'loopcnt %06d carr %+010.2f  pd%.3f\n',...
                    loopcnt,carrier_dco_deltaL*Freq_sample-Fc,pd);
                fclose(fid);
                %                 fprintf('P_k对角线: [%.2e, %.2e, %.2e, %.2e]\n', diag(P_kf));                                 
            end
            %% 本地载波更新已完成
            %% 计算本地与接受多普勒差值
            FCARR(loopcnt)=carrier_dco_deltaL*Freq_sample-Fc;
            freceive=yv_new(mstime);
            doperror=freceive-carrier_dco_deltaL*Freq_sample+Fc;
            DOPERRORre(loopcnt)=doperror;
            %% 判断失锁？
            if abs(doperror) > doppler_fail_threshold
                bad_ms = bad_ms + 1;
            else
                bad_ms = 0;
            end
            if bad_ms >= consecutive_ms_for_fail
                fprintf('跟踪失败：|Δf|>%g Hz 连续 %d ms（mstime=%d）。\n', ...
                    doppler_fail_threshold, consecutive_ms_for_fail, mstime);
                success = 0;
                break;  % 只跳出 while，本次 Monte-Carlo 结束
            end
            if(mode~=5)
                %% 码环路
                if( mod(mscountL1_bitperiod,codeloop_cohtime) == 0)
                    codeloop_noncohcnt = codeloop_noncohcnt + 1;
                    corrnoncohsum_code = corrnoncohsum_code + abs(corrcohsum_code).^2;
                    corrcohsum_code =  zeros(1,LengthCorreArray);
                    if(codeloop_noncohcnt == codeloop_noncohtime)
                        codeloop_noncohcnt = 0;
                        E = corrnoncohsum_code(mid - 5);
                        L = corrnoncohsum_code(mid + 5);
                        cd = (E-L)/(E+L);
                        vo = vo_1 + cd1*(cd - cd_1) + cd2*(cd + cd_1);  %VCO输入
                        code_dco_deltaL = code_dco_deltaL - (vo - vo_1)/Freq_sample; %本地复现码频率更新
                        %                 code_dco_deltaL*Freq_sample-Freq_code
                        %本地复现码频率存入数组
                        
                        vo_1 = vo;
                        cd_1 = cd;
                        corrnoncohsum_code =  zeros(1,LengthCorreArray);
                    end
                end
                %% 载噪比估计
                if( mod(mstime,jugujitime) == 0)
                    %该载噪比估计不适用！
                    
                    % %             if(cntL >= 2)
                    %                 tmmp = accum;
                    %                 NBP = abs(sum(tmmp))^2;
                    %                 I2_Q2 = real(sum(tmmp))^2-imag(sum(tmmp))^2;
                    %                 WBP = tmmp * tmmp';
                    %                 NP = NBP/WBP;
                    %                 loopcnt;
                    %                 tt1=10*log10((NP-1)/(20-NP)*1000);
                    %                 pld = I2_Q2/NBP;
                    % %                 carrier_dco_deltaL*Freq_sample-Fc
                    % %                 code_dco_deltaL*Freq_sample-Freq_code
                    %                 NParray = [NParray NP];
                    %                 PLD = [PLD pld];
                    %                 PD = [PD pd];
                    %                 accum = [];
                    %                 CN0L=[CN0L tt1];
                    %                 cntL = 0;
                    %
                    %                 fid = fopen(filename,'a');
                    %                 fprintf(fid,'loopcnt %06d carr %+010.2f code %+010.2f NP %+010.2f PLD %+010.2f\n',...
                    %                     loopcnt,carrier_dco_deltaL*Freq_sample-Fc,code_dco_deltaL*Freq_sample-Freq_code,NP,pld);
                    %                 fclose(fid);
                    % %             end
                    
                    %改成矩估计
                    Tju=jugujitime/2;
                    M2CO=M2/Tju;
                    M4CO=M4/Tju;
                    tt1=10*log10((sqrt(2*M2CO^2-M4CO))/(M2CO-sqrt(2*M2CO^2-M4CO))*1000);
                    M2=0;
                    M4=0;
                    
                end
            end
            %% 本地码 CNO估计完成
            FCODE(loopcnt) =  code_dco_deltaL*Freq_sample-Freq_code;
            CN0L(loopcnt)=tt1;
            
        end
        %%
        %信号生成模块
        result = randn*sigma + sqrt(-1)*randn*sigma; %%%%%%%% 生成IQ两路信号噪声
        %%%%%%%% 低中频到零中频变换得到最终结果 %高速电文以及映射至CA码里。
        for ss= 1:4
            result = result + A(ss)*codelow(ss,chipcounter(ss)+1)*(1-gt(ss))*BITlow(bitnumlow(ss)+1)*exp(sqrt(-1)*2*pi*carrier_dco_phase(ss))+...
                A(ss)*codehigh(ss,chipcounter(ss)+1)*gt(ss)*exp(sqrt(-1)*2*pi*carrier_dco_phase(ss));
            
        end
        tmptmp = floor(mod(chipcntL+code_dco_phaseL+LocalCorreArray+codelength,codelength)+1); %本地扩频码组对应的码片索引
        tmp1 = LocalCode(tmptmp);   %由本地扩频码组对应的码片索引，得到对应码片
        
        csktmptmp=floor(mod(chipcntL+code_dco_phaseL+codelength,codelength)+1);
        tmpcskcode=Localhightraverse(:,csktmptmp);
        %     corrsum = corrsum + result*exp(-sqrt(-1)*2*pi*carrier_dco_phase(AcqSatNUM))*tmp1;
        if(mod(mscount(AcqSatNUM),2)==0)
            if(mode==3||mode==6||mode==9)
                corrsum = corrsum + result*exp(-sqrt(-1)*2*pi*carrier_dco_phaseL)*tmp1*BITlow(bitnumlow(AcqSatNUM)+1);%mode==3 EKF分支，IQ观测值需要剥离电文
            elseif(mode==5)%完全对齐
                tmptmp = floor(mod(chipcounter(AcqSatNUM)+code_dco_phase(AcqSatNUM)+LocalCorreArray+codelength,codelength)+1);
                tmp1 = LocalCode(tmptmp);
                csktmptmp=floor(mod(chipcounter(AcqSatNUM)+code_dco_phase(AcqSatNUM)+codelength,codelength)+1);
                tmpcskcode=Localhightraverse(:,csktmptmp);
                corrsum = corrsum + result*exp(-sqrt(-1)*2*pi*carrier_dco_phase(AcqSatNUM))*tmp1;
            else
                corrsum = corrsum + result*exp(-sqrt(-1)*2*pi*carrier_dco_phaseL)*tmp1;
            end
        else
            %完成CSK解调的相关 64 组相关值
            if(demodulation)
                corrsumforcskdemodulate=corrsumforcskdemodulate+result*exp(-sqrt(-1)*2*pi*carrier_dco_phaseL).*tmpcskcode';
            end
        end
        n = n + 1;
        
        for i = 1:4
            carrier_dco_phase(i) = carrier_dco_phase(i) + carrier_dco_delta(i);
            code_dco_phase(i) = code_dco_phase(i)+ code_dco_delta(i) + 0;
            
            if(carrier_dco_phase(i) > 1)
                carrier_dco_phase(i) = carrier_dco_phase(i) - 1;
            end
            if(code_dco_phase(i) > 1)
                code_dco_phase(i) = code_dco_phase(i) - 1;
                chipcounter(i) = chipcounter(i) + 1;
                %%% 伪码长度codelength个码片，若计数器计到，则经过一个伪码周期，计数器归零，加一毫秒
                if(chipcounter(i) >= codelength)
                    chipcounter(i) = 0;
                    mscount(i) = mscount(i) + 1;
                    %% 每一毫秒改变DOP 重新计算多普勒/码步进
                    index=ceil(n/(Freq_sample*0.001))+1;
                    carrier_dco_delta(i)=(Fc+yv_new(index))/Freq_sample;
                    code_dco_delta(i)=(yv_new(index)/codedopcoef+Freq_code)/Freq_sample ;
                    if(mscount(i) == 4)
                        bitnumlow(i) = bitnumlow(i) + 1;
                        mscount(i) = 0;
                        gt(i) = 0;
                    elseif(mscount(i) == 1)
                        xk(i) = 32*BIThigh(bitnumhigh(i)+1)+16*BIThigh(bitnumhigh(i)+2)+8*BIThigh(bitnumhigh(i)+3)+4*BIThigh(bitnumhigh(i)+4)+...
                            2*BIThigh(bitnumhigh(i)+5)+BIThigh(bitnumhigh(i)+6);%注意BITHIGH值为0或1
                        codehigh(i,:) = circshift(codehigh_initial(i,:),xk(i));  %检查CA码是否移位！！！
                        bitnumhigh(i) = bitnumhigh(i) + 6;
                        gt(i) = 1;
                    elseif(mscount(i) == 2)
                        %低速每个符号发两个低速扩频码周期，bitnumlow不需要+1
                        gt(i) = 0;
                    elseif(mscount(i) == 3)
                        xk(i) = 32*BIThigh(bitnumhigh(i)+1)+16*BIThigh(bitnumhigh(i)+2)+8*BIThigh(bitnumhigh(i)+3)+4*BIThigh(bitnumhigh(i)+4)+...
                            2*BIThigh(bitnumhigh(i)+5)+BIThigh(bitnumhigh(i)+6);%注意BITHIGH值为0或1
                        codehigh(i,:) = circshift(codehigh_initial(i,:),xk(i));  %检查CA码是否移位！！！
                        bitnumhigh(i) = bitnumhigh(i) + 6;
                        gt(i) = 1;
                    else
                        disp('err');
                        pause;
                    end
                end
                
            end
        end
        
        
        %% 本地复现信号的调整如下
        localn = localn + 1;
        carrier_dco_phaseL = carrier_dco_phaseL + carrier_dco_deltaL ;
        if(carrier_dco_phaseL > 1)
            carrier_dco_phaseL = carrier_dco_phaseL - 1;
        end
        code_dco_phaseL = code_dco_phaseL + code_dco_deltaL;
        
        if(code_dco_phaseL > 1)
            code_dco_phaseL = code_dco_phaseL - 1;
            chipcntL = chipcntL + 1;
            if(chipcntL >= codelength)
                dump = 1;
                chipcntL = 0;
                mscountL = mscountL + 1;
                mstime=mstime+1;
                %%% 电文速率50bps，若计数器计到，则经过一个周期，对应一个电文码，毫秒计数器归零，比特计数器加一
                if(mscountL >= bitperiod)
                    mscountL = 0;
                    %                 bitnumL = bitnumL+1;
                    
                end
            end
        end
        
        if(mod(localn,Freq_sample * 0.001) == 0)          %每1ms测量一下跟踪误差
            %         localn = 0;
            cod = code_dco_phaseL - code_dco_phase(AcqSatNUM);
            car = carrier_dco_phaseL - carrier_dco_phase(AcqSatNUM);
            codeerr(loopcnt+1) =  cod;
            carrerr(loopcnt+1) =  car;
        end
        
        if(n > Freq_sample*runtime)
            %        figure, stem(histnumi);
            %        figure, stem(histnumq);
            break
        end
        
        %%%%%%%% sec加一秒并打印
        samplenum = samplenum + 1;
        if(samplenum == Freq_sample)
            samplenum = 0;
            sec = sec + 1;
            %         fprintf('%ds....\n',sec);
            timeold=t;
            t = toc;  % 获取当前耗时
            thissecond=t-timeold;
            waitbar(sec / runtime, h, sprintf('第%d/%d秒数据 本帧耗时%.f秒 总耗时%.f秒', sec, runtime,thissecond, t));
        end
    end
    if exist('h','var') && ishghandle(h), close(h); end  % 稳妥关闭进度条
    runtime_s = toc(t0);                                  % 本次用时(秒)
    
    % —— 失败：记录并进入下一次 MC ——
    if ~success
        mc_log{mc, :} = [seeds(mc), CN0, mode, mstime, ...        % seed, CN0, mode, fail_ms
            NaN, NaN, NaN, NaN, NaN, NaN, ...        % BER/RMSE/CN0mean占位
            runtime_s, 0];                           % 用时, success=0
        continue;   % 直接进入下一次 Monte-Carlo
    end
    if(success==1)
        disp("跟踪结束，结果统计：");
        %% 统计BER
        if(demodulation)
            highbitcskcosysafteerLDPC = reshape(decode_matrix(:, 1:decode_index-1), 1, []);
            [BER_coherent,BER_non_coherent,BER_LDPC_coherent]=countBERandSER(runtime,BITHIGHTRANS,highbitcskco_symbol,highbitcsknoco_symbol,highbitcskcosysafteerLDPC);
        else %设定误码率无效值
            BER_coherent=-10;
            BER_non_coherent=-10;
            BER_LDPC_coherent=-10;
        end
        %% 统计跟踪误差
        if(mode~=5)
            [carrier_phase_RMSE,doppler_RMSE,CN0MEAN]=Trackplotandcount(runtime,carrerr,DOPERRORre,CN0L,FCODE,FCARR,codeerr);
        end
        mc_log{mc, :} = [seeds(mc), CN0, mode, NaN, ...                 % 无 fail_ms
            BER_coherent, BER_non_coherent, BER_LDPC_coherent, ...
            carrier_phase_RMSE, doppler_RMSE, CN0MEAN, runtime_s, 1];
        
        % 获取工作区中的所有变量
        workspaceVariables = whos;  % 获取所有变量的详细信息
        % 生成文件名，使用两个整数变量作为文件名的一部分
        fileName = sprintf('workspace_data_%d_%d_time%d.mat', CN0, mode,mc);
        % 将所有工作区的变量保存到 .mat 文件中
        save(fileName);
        disp(['工作区的所有数据已成功保存为 .mat 文件：' fileName]);
    end
    % fclose(fid);
end
writetable(mc_log, 'simulation_results_MC.xlsx', 'Sheet','MC','WriteMode','overwritesheet');
save('MC_all_runs.mat','mc_log');
disp('Monte-Carlo 全部完成。');