% TCD 浮式风力涡轮机代码
% 由 Saptarshi Sarkar 在 2017-2019 年间作为其博士项目开发
% 目的：模拟 5MW OC3 浮式风力涡轮机的动态行为，包括空气动力学、结构动力学和系泊系统交互。

clc; clear -globals; clear;
clear SystemMatrices BaselineControllers;
% 清除命令窗口、工作区和全局变量，确保运行环境干净。
% 清除预定义的系统矩阵和基线控制器变量。

addpath('./Wind_Dataset','./Wave_Dataset','./src') 
% 将包含风场数据、波浪数据和源代码的目录添加到 MATLAB 搜索路径。

pause(.5);
% 暂停执行 0.5 秒，以确保库加载或初始化过程顺利完成。

global GRAVACC FLUIDDENSITY;      
GRAVACC = 9.80655; FLUIDDENSITY = 1025;
% 定义全局常量：
% - GRAVACC：重力加速度 (m/s^2)。
% - FLUIDDENSITY：海水密度 (kg/m^3)，用于波浪和系泊系统计算。

% 卸载之前加载的系泊动力学库，以避免冲突。
if libisloaded('MoorDyn') || libisloaded('MoorApiwin64')
    if libisloaded('MoorApiwin64')
        calllib('MoorApiwin64','finish'); 
        unloadlibrary('MoorApiwin64');
    else
        calllib('MoorDyn','LinesClose');      
        unloadlibrary MoorDyn;  
    end                                    
    % 确保正确清理 MoorDyn 或 MoorApiwin64 动态链接库。
end

% 读取风力涡轮机模型的输入数据文件。
[Airfoils, Geometry] = ReadWindTurbineAeroData('rad');
% 读取空气动力学数据（翼型特性和叶片几何），角度单位为弧度。
[Blade, Tower]       = ReadWindTurbineStructuralData();
% 读取叶片和塔架的结构特性数据。
[ElastoDyn]          = ReadElastoDyn();
% 读取 ElastoDyn 输入文件，包含风力涡轮机的运动学和初始条件。
[Servo]              = ServoDyn();
% 读取 ServoDyn 输入文件，包含控制系统参数（如桨距和扭矩控制）。

% 创建塔架和叶片结构，用于仿真。
Twr = CreateTwr(Tower, ElastoDyn); 
% 根据塔架数据和 ElastoDyn 参数创建塔架结构。
Bld = CreateBld(ElastoDyn, Geometry, Blade); 
% 根据 ElastoDyn、几何和叶片数据创建叶片结构。
Platform = CreatePlatform();
% 创建浮式平台结构。
Platform.Mooring     = 1;    % 1 表示使用 MoorDyn 系泊模型，2 表示使用 OpenMoor。
Platform.WaveLoads   = 0;    % 1 表示使用 Morison 方程计算浮筒上的波浪载荷，0 表示不计算。

% 定义可用的自由度（DOFs）。
DOFsStr = {'Sg','Sw','Hv','R','P','Y','TFA1','TSS1','TFA2','TSS2','NacYaw','GeAz','DrTr','B1F1','B1E1','B1F2','B2F1','B2E1','B2F2','B3F1','B3E1','B3F2'};
% 自由度字符串列表，表示平台、塔架和叶片的运动模式（如平移、旋转、塔架前倾、叶片挥舞等）。

% 处理波浪载荷选项。
if Platform.WaveLoads == 1
    waveopt.wave_profile = 'load';
    waveopt.wave_file = 'WaveHs3_Tp11_Dir90_NoCur0_Tf50';
    wave = Wave(waveopt);   
    % 如果启用波浪载荷，加载指定的波浪数据文件（Hs=3m, Tp=11s, 波向90度，无洋流，模拟时间50s）。
else
    wave.wave_file = [];
    % 如果不计算波浪载荷，设置波浪文件为空。
end

% 定义关闭的自由度。
TurnedOffDOFsStr = {};
% 列出需要关闭的自由度（此处为空，表示所有自由度都启用）。
TurnedOffDOFs   = find(ismember(DOFsStr, TurnedOffDOFsStr));
% 查找关闭自由度的索引。
DOFs.Avail      = length(DOFsStr);
% 可用自由度总数。
DOFs.Active     = setdiff(1:length(DOFsStr), TurnedOffDOFs);
% 激活的自由度（排除关闭的自由度）。
DOFs.ActNominal = setdiff(1:22, TurnedOffDOFs);
% 标称激活自由度（前22个自由度中排除关闭的）。
DOFs.nDOFs      = length(DOFs.Active);
% 激活自由度的数量。

nDOFs = length(DOFs.Active);
% 存储激活自由度的总数。
if find(ismember(TurnedOffDOFsStr, 'GeAz'))
    q0 = zeros(1, 2*nDOFs);
    % 如果齿轮箱转角（GeAz）关闭，初始状态向量全设为零。
else
    q0 = zeros(1, 2*nDOFs);
    q0((12 - sum(TurnedOffDOFs < 12)) + nDOFs) = ElastoDyn.RotSpeed;
    % 否则，设置初始状态向量，其中转子转速（GeAz）设为 ElastoDyn 中定义的初始转速。
end

% 读取 TurbSim 生成的风场数据和网格。
[velocity, Wind.y, Wind.z, Wind.nz, Wind.ny, Wind.dz,...
          Wind.dy, Wind.dt, Wind.zHub, Wind.z1, Wind.SummVars] = readBLgrid('NTM_A');
% 读取 TurbSim 风场文件（NTM_A），获取速度场、网格点坐标、时间步长等。
Wind.t_TurbSim = (0:1:size(velocity,1)-1)*Wind.dt;
% 计算风场的时间向量。

% 如果需要使用稳态风场，取消注释以下行。
% [Wind.t_TurbSim, velocity] = SteadyWind(Wind.y, Wind.z, 11.4, 0.2);
% 生成均匀风速为 11.4 m/s，湍流强度为 0.2 的稳态风场。

% 创建风速插值器，用于动态计算任意点和时间处的风速。
gv = {Wind.t_TurbSim, [1 2 3], Wind.y, Wind.z};
Wind.Velocity = griddedInterpolant(gv, velocity, 'linear');
% 使用线性插值方法生成风速插值器。

% 设置风场和控制参数。
WindNom.PittandPeters = true(0);  % 是否启用 Pitt and Peters 模型（0 为禁用）。
WindNom.PitchControl  = true(1);  % 是否启用桨距控制（1 为启用）。
WindNom.AeroElastic   = true(1);  % 是否考虑顺桨叶方向的空气弹性效应（1 为启用）。
WindNom.y = Wind.y;
WindNom.z = Wind.z;
% 存储风场网格的 y 和 z 坐标。

% 初始化流入角。
phi = InitializeInflowAngle(ElastoDyn, Bld, Wind);
% 根据 ElastoDyn、叶片和风场数据初始化叶片的流入角。

% 初始化系泊系统。
mooring_load_ptr = zeros(1,6);           % 创建系泊力指针，初始为零（6个自由度）。
mooring_load_ptr = libpointer('doublePtr', mooring_load_ptr);  % 将系泊力数组转换为指针，供外部函数修改。

if Platform.Mooring == 1
    [notfound, warnings] = loadlibrary('MoorDyn', 'MoorDyn');     % 加载 MoorDyn 动态链接库。
    mooring_status = calllib('MoorDyn', 'LinesInit', zeros(6,1), zeros(6,1));   % 初始化 MoorDyn 系泊系统。
    if mooring_status ~= 0
        error('MoorDyn 无法加载');
    else
        disp('MoorDyn 加载成功');
    end
elseif Platform.Mooring == 2
    loadlibrary('MoorApiwin64', 'moorapi.h'); 
    input_file = 'CaseOC3.xml'; 
    calllib('MoorApiwin64', 'initialize', input_file);
    % 加载 OpenMoor 库并使用 OC3 配置文件初始化。
else
    error('请选择正确的系泊动力学模型');
end

% 初始化控制变量。
Controls0 = [Servo.VS_RtTq, ElastoDyn.BlPitch, reshape(phi, [1 3*Bld.nb])];
% 初始控制变量包括：额定扭矩（Servo.VS_RtTq）、初始桨距角（ElastoDyn.BlPitch）和展向流入角。

% 设置仿真时间参数。
t0  = 0;      % 仿真起始时间 (s)。
tf  = 5;      % 仿真终止时间 (s)。
deltat = 0.0125;  % 时间步长 (s)。
t = t0:deltat:tf; % 时间向量。
n = (tf-t0)/deltat; % 时间步数。

% 运行 ode4（四阶定步长 Runge-Kutta 方法）进行动态仿真。
[q, Controls] = ode4(@(t,q,Controls) RHS(t, q, Controls, DOFs, ElastoDyn, Airfoils, Twr, Bld,...
                Platform, Wind, WindNom, wave, mooring_load_ptr, Servo), t0, deltat, tf, q0, Controls0);
% 调用 RHS 函数计算状态导数，q 为状态向量，Controls 为控制向量。

% 清理系泊系统库。
if Platform.Mooring == 1
    calllib('MoorDyn', 'LinesClose');       % 关闭 MoorDyn。
    unloadlibrary MoorDyn;                  % 卸载 MoorDyn 库。
elseif Platform.Mooring == 2
    pause(1)
    calllib('MoorApiwin64', 'finish');
    unloadlibrary('MoorApiwin64');
    % 关闭并卸载 OpenMoor 库。
end

% 从仿真状态生成结果。
[OopDefl1, IPDefl1, TipDzb1, TTDspFA, TTDspSS, LSSTipV, NacYaw, TwrClrnB1, TwrClrnB2,...
                    TwrClrnB3, BldTwrHrClr, Surge, Sway, Heave, Roll, Pitch, Yaw] = Results(t, q, Controls, DOFs, ElastoDyn, Twr, Bld);
% 计算并输出仿真结果，包括：
% - OopDefl1, IPDefl1: 叶片挥舞和顺桨方向的变形。
% - TipDzb1: 叶片尖端位移。
% - TTDspFA, TTDspSS: 塔架前倾！

%% Plots
figure(1)
hold on
plot(t,OopDefl1)
xlabel('Time')
ylabel('OopDefl1')

figure(2)
hold on
plot(t,IPDefl1)
xlabel('Time')
ylabel('IPDefl1')


figure(3)
hold on
plot(t,LSSTipV)
xlabel('Time')
ylabel('LSSTipV')

figure(4)
hold on
plot(t,NacYaw)
xlabel('Time')
ylabel('NacYaw')
 
figure(5)
hold on
plot(t,TTDspFA)
xlabel('Time')
ylabel('TTDspFA')
 
figure(6)
hold on
plot(t,TTDspSS)
xlabel('Time (sec)','FontSize',12)
ylabel('Tower side-to-side displacement (m)','FontSize',12)
 
figure(7)
hold on
plot(t, Surge)
xlabel('Time')
ylabel('Platform Surge')

figure(8)
hold on
plot(t, Sway)
xlabel('Time')
ylabel('Platform Sway')

figure(9)
hold on
plot(t, Heave)
xlabel('Time')
ylabel('Platform Heave')

figure(10)
hold on
plot(t, Roll)
xlabel('Time')
ylabel('Platform Roll')

figure(11)
hold on
plot(t, Pitch)
xlabel('Time')
ylabel('Platform Pitch')

figure(12)
hold on
plot(t, Yaw)
xlabel('Time')
ylabel('Platform Yaw')

figure(13)
plot(t, TwrClrnB1, t, TwrClrnB2, t, TwrClrnB3)
xlabel('Time')
ylabel('Tower Blade Clearence')

figure(14)
hold on
plot(t, Controls(:,2), t, Controls(:,3), t, Controls(:,4))
xlabel('Time')
ylabel('Blade pitch angles')
