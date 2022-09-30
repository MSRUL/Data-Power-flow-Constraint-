clear,clc
close all;
path_base = mfilename('fullpath');i=findstr(path_base,'\');path_base=path_base(1:i(end));
path_last = mfilename('fullpath');i=findstr(path_base,'\');path_last=path_last(1:i(end-1)); %此代码必须运行m文件执行
%% ZYX 202203


%% 【！！】Pay attention to the parameters before modification
%------------Case Name------------------------------
Case_Name='IEEE141'; %'IEEE33' 'IEEE69' 'IEEE85' 'IEEE141'
Opt_Variable='Reactive_Power'; %'Active_Power' 'Reactive_Power'
MN_info='MN_no';   % MN_no   ME_in  %  Measure_Noise
IDL_Peremeters.type='IDL_QuaModified';  % 'IDL'  'IDL_QuaModified'
Objective_Info.type='Min Adjustment';  %'Min Loss'  'Min Adjustment' 'Hybrid' 'Min_Voltage_deviation_rate' 'Min_VDR_and_Min_Adj' 'DG_Consumption'
%%
start_in_1000=571;


switch Case_Name
    case 'IEEE33'   
        case_name_or=case33bw;
    case 'IEEE69'
        case_name_or=case69;
    case 'IEEE85'   
        case_name_or=case85;
    case 'IEEE141'
        case_name_or=case141;
end

%------------Generated Data Parameters--------------
Range=3.2;   % 3.2 In Letter Cases
New_Training_data=0;
% Only When New_Training_data==1:
Measure_Noise.range=0;  %Only New_Training_data==1
training_data_file='1000组_Range3.2_Vm_Noise0_no1'; %Only New_Training_data==0  2_Vm_Noise0_Col13833_no1  2_Vm_Noise0_no1
testing_data_file='1000组_Range3.2_Vm_Noise0_no2'; %Only New_Training_data==0
% 1000组_Range3.2_Vm_Noise0_no1
% 1000组_Range3.3_Vm_Noise0_no1

%------------Dimension Lifting Parameters-----------
num_inpuit=1000;   %Num of Training Data Sets--------
Nrbf1=10;  % 10 In Letter Cases
Nrbf2=0; %input_num-output_num;
rbf_type='polyharmonic';rbf_type_Inverse='polyharmonic';
% 'polyharmonic'  'invmultquad' 'invquad'  '3n+1'  '2second order'
% 'MIX polyh invquad invmultquad'  'MIX polyh invquad invmultquad 3n+1'
%'MIX invquad and invmultquad'  'MIX polyh and invmultquad' 'MIX polyh and invquad'
%'MIX polyh and 3n+1'
% 采用polyh invquad invmultquad三种混合的升维时，设定的升高维度自动升高1.5倍

%%
%DD_PF Mode【DDPF Type】
% IDL_Peremeters.type='IDL';  % 'IDL'  'IDL_withSM'  'IDL_QuaModified'  'IDL_QuaModified_withSM'  'IDL_QuaModified_withSM_3Np1'  'IDLSM2'  'IDLSM3'  'IDL_CubicSpace'
%IDL                            :1: incomplete dimension lifting;
%IDL_withSM                     :2: The incomplete dimensiong lifting with senstive variable
%IDL_QuaModified                :3: the basic incomplete dimension lifting with the control variable is quadratic space（2N+1）
%IDL_QuaModified_withSM         :4: The incomplete dimensiong lifting with senstive variable with the control variable is quadratic space（2N+1）
%IDL_QuaModified_withSM_3Np1    :5: The incomplete dimensiong lifting with senstive variable with the control variable is quadratic space（Entire 3N+1）
%IDLSM2                         :6: The simplified IDL_withSM    ；Senstive Matrix of u is the linear function of the lifted x        
%IDLSM3                         :7: The simplified IDLSM2        ；Senstive Matrix of u is the linear function of the Unlifted x  (orginal x)  
%IDL_CubicSpace                 :8: The Cubic Space 3Np1       

% --------Output Variable List----------------------
ItemName='';
TestItem.PQ_Vm=1;       if TestItem.PQ_Vm ;ItemName=strcat(ItemName,'Vm'); end
TestItem.Va=0;          if TestItem.Va ;ItemName=strcat(ItemName,'Va'); end
TestItem.LP=0;          if TestItem.LP ;ItemName=strcat(ItemName,'LP'); end
TestItem.LQ=0;          if TestItem.LQ ;ItemName=strcat(ItemName,'LQ'); end
TestItem.PLoss=0;       if TestItem.PLoss ;ItemName=strcat(ItemName,'PLoss'); end

if New_Training_data
    KPM_Parameters_Info=strcat(num2str(num_inpuit),'组_Range',num2str(Range),'_',ItemName,'_Noise',num2str(Measure_Noise.range));
else
    KPM_Parameters_Info=training_data_file;
end
additional_info_KMP='Using_M_equa_0__';
additional_info_KMP=strcat(additional_info_KMP,rbf_type);
KPM_Parameters_Info=strcat(KPM_Parameters_Info,additional_info_KMP);


%------------Opt_Variable---------------------------

% --------Objective in Optimization List----------------------
% see details in process below
%% Record the case info (Opt info) The deteail of the parameters are not included.
case_obj_info=Case_Name;
case_obj_info=strcat(case_obj_info,'_',IDL_Peremeters.type);
switch Objective_Info.type
    case 'Min Adjustment'
        case_obj_info=strcat(case_obj_info,'_MinAdj');
    case 'Min_Voltage_deviation_rate'
        case_obj_info=strcat(case_obj_info,'_MinVDR');
end
additional_info_IDL='';
case_obj_info=strcat(case_obj_info,additional_info_IDL);

switch Opt_Variable
    case 'Reactive_Power'
        ;
    case 'Active_Power'
        case_obj_info=strcat(case_obj_info,'_OptP');
end
switch MN_info
    case 'MN_no'
        ;
    case 'ME_in'
%        Sigma=0.0005;  % About 95% of the data is less than 0.1% error
%         case_obj_info=strcat(case_obj_info,'_Noised');
       Sigma=0.005;  % About 95% of the data is less than 1% error
        case_obj_info=strcat(case_obj_info,'_Noised2'); 
end

if strcmp(IDL_Peremeters.type,'IDL_QuaModified')
    Additional_Info='QS_Test_01';
else
    Additional_Info='';
end

if start_in_1000~=1
Additional_Info=strcat(Additional_Info,'_',num2str(start_in_1000),'start');
end
case_obj_info=strcat(case_obj_info,Additional_Info); 

%% Init the Case Parameters
case_name_or_simp=Transfer_node_num_to_consecutive(case_name_or);  % Change the bus name to consecutive numbers
case_name=Transfer_node_num_to_consecutive(case_name_or);  % Change the bus name to consecutive numbers
Bus_Num=size(case_name.bus,1);
line_active=find(case_name.branch(:,11)==1);
Branch_Num=size(line_active,1);
Standard_result=runpf(case_name);
fprintf('OR Power Loss: %f\n',sum(abs(Standard_result.branch(:,14)+Standard_result.branch(:,15)))/Standard_result.baseMVA);
[ref,pv, pq,pv_pos,ref_pos,pv_total] = Case_Info(case_name);
% figure;plot(Standard_result.bus(:,8))
%% Init the PhotoVoltak/CapacitorBank Parameters
HaveDG=1; %This Program default is 1(For the optimization)
switch Case_Name
    case 'IEEE33'
        [Device_Info] = Init_Device_Info(case_name,HaveDG);
    case 'IEEE69'
        [Device_Info] = Init_Device_Info_Case69(case_name,HaveDG);
    case 'IEEE85'
        [Device_Info] = Init_Device_Info_Case85(case_name,1);
    case 'IEEE141'
        [Device_Info] = Init_Device_Info_Case141(case_name,HaveDG);
end
%% Data Parameter
input_num=Device_Info.Input_Num;
test_input_num=input_num; %size(pq,1)*2+size(pv,1);
case_name.bus(pq,3:4)=0;
case_name.gen([pv_total;ref_pos],[2,6])=0;
num_test_input=max(num_inpuit*1.2,num_inpuit+1000);
% num_test_input=1;
%% Position of Variable in Input Vector
[Pos_In_OutPut] = Define_Pos_VariableVector(case_name,Device_Info,Bus_Num,Branch_Num);
%--------Additional Modification---------------------
switch Opt_Variable
    %  Pos_In_OutPut.Pos_Q_in_Inventer  Pos_In_OutPut.Pos_P_in_Inventer  [Pos_In_OutPut.Pos_P_in_Inventer;Pos_In_OutPut.Pos_Q_in_Inventer]
    case 'Reactive_Power'
        fprintf('\n无功功率优化');
        Pos_In_OutPut.Un_RisePos=[Pos_In_OutPut.Pos_Q_in_Inventer];% Position of the Unrised Vector 【Control Variable】
    case 'Active_Power'
        fprintf('\n有功功率优化');
        Pos_In_OutPut.Un_RisePos=[Pos_In_OutPut.Pos_P_in_Inventer];% Position of the Unrised Vector 【Control Variable】
    case 'DG_PG_Power'
        Pos_In_OutPut.Un_RisePos=[Pos_In_OutPut.Pos_P_in_Inventer;Pos_In_OutPut.Pos_Q_in_Inventer];% Position of the Unrised Vector 【Control Variable】
end
Pos_In_OutPut.RisePos=Find_Wihtout_Which([1:input_num],Pos_In_OutPut.Un_RisePos); % Position of the Rised Vector 【State Variable】
Pos_In_OutPut.Pos_Input_Rise=length(Pos_In_OutPut.RisePos)+1:Pos_In_OutPut.input_num; %The Unreised Part of Input During the Process of maltiping the M.

%% Generate the Training Data
Measure_Noise.state.input=1;% Use the noised data
Measure_Noise.state.output=1;% Use the noised data
[Tranning_Range,Testing_Range] = Init_Training_Data_Range(Range);
% ------------------------------------------
% Nrbf1=100;testnum=10;
% RecordRank=TestRank(testnum,input_num,num_inpuit,Tranning_Range,Device_Info,case_name,case_name_or_simp,'polyharmonic', rand(test_input_num,Nrbf1)); % 'MIX polyh and 3n+1'
% ------------------------------------------
if New_Training_data
    %     [Input,Output] = Generate_Tranning_OR_Testing_Data(input_num,test_input_num,num_inpuit,Tranning_Range,Device_Info,case_name,case_name_or_simp,Measure_Noise,TestItem);
    %%  Storge Data
    for i=1:5
        [Input,Output] = Generate_Tranning_OR_Testing_Data(input_num,test_input_num,num_inpuit,Tranning_Range,Device_Info,case_name,case_name_or_simp,Measure_Noise,TestItem);
        Data.Input=Input;
        Data.Output=Output;
        save_data_path=strcat(path_base,'【变量】\',Case_Name,'训练数据\',num2str(num_inpuit),'组_Range',num2str(Range),'_',ItemName,'_Noise',num2str(Measure_Noise.range),'_no',num2str(i),'.mat');
        save(save_data_path,'Data');
    end
    training_data_file=strcat(path_base,'【变量】\',Case_Name,'训练数据\',num2str(num_inpuit),'组_Range',num2str(Range),'_',ItemName,'_Noise',num2str(Measure_Noise.range),'_no',num2str(1),'.mat');
else
    load(strcat(path_base,'【变量】\',Case_Name,'训练数据\',training_data_file,'.mat'));
%     Data=Data_CoL; % When collinear data is in
end
% ------------------------------------------

if strcmp(MN_info,'ME_in')  % The Measure Error is considered in the data Training stage
     Data.Input= Data.Input.*(normrnd(1,Sigma,size(Data.Input)));
     Data.Output= Data.Output.*(normrnd(1,Sigma,size(Data.Output)));
end

Input_Test_vif=Data.Input;
pos=find(Input_Test_vif(:,1)==0);
Input_Test_vif(pos,:)=[];
pos=find(Input_Test_vif(:,1)==1);
Input_Test_vif(pos,:)=[];
[VIF] = VIF_Cal(Input_Test_vif);
fprintf('\n VIF:%f\n',VIF);
%%
% Data.Input=Data.Input(:,1:120);
% Data.Output=Data.Output(:,1:120);
%% Koopman Matrix Calculate
fprintf('\n\nUse the Training Data to Test:\n');

% [M_Compelet,M_Inverse]=Cal_M_Koopman(Input,Output,rbf_type,rbf_type_Inverse,cent_input_compelete,cent_output);
%-----------【Incompelet】---------------------------------   The【Compelet】can be find in Koopman_PF_Optimization.m
IDL_Peremeters.rbf_type=rbf_type;
IDL_Peremeters.RisePos=Pos_In_OutPut.RisePos;
IDL_Peremeters.Un_RisePos=Pos_In_OutPut.Un_RisePos;
IDL_Peremeters.cent = Generate_Cent(IDL_Peremeters,Nrbf1);
IDL_Peremeters.type;
IDL_Peremeters.M_Asstence=rand(1,length(Pos_In_OutPut.Un_RisePos));

 Train_Ops.way='LS_OR';   %  'LS_OR'  'LS_OR_pu'  'LS_Opt' 'Load_Trained_Data'
 Train_Ops.Beita=0.1;
 Train_Ops.M_Type=1; % Default
if strcmp(Train_Ops.way,'Load_Trained_Data')
    % So far ; 【Only】 the 3Np1 type is meaningful. 
    load(strcat(path_base,'【变量】\',Case_Name,'训练M矩阵结果\','M_Opt_withSens__Alpf_is_0.1.mat'));% M
    M_Incompelet=M.M_Incompelet;M_Struct=M.M_Struct;

    [Output_Test_KPMCaled] = Case_PF_Accuracy_Test(IDL_Peremeters,M_Incompelet,M_Struct,Data.Input);
    fprintf('\nTest DD IDL PF， Average Error： %f   ;Maximum Error:  %f  ',mean(mean(abs(Output_Test_KPMCaled-Data.Output))),max(max(abs(Output_Test_KPMCaled-Data.Output))));

else
    %% The Training Parts
   switch IDL_Peremeters.type % 'IDL' 'IDL_QuaModified'
        case 'IDL'
            [M_Incompelet,M_Struct,M_Output] = Cal_KMP_Total(Data,IDL_Peremeters,Train_Ops);
        case 'IDL_QuaModified'
            [M_Incompelet,M_Struct,M_Output] = Cal_KMP_Total_QS(Data,IDL_Peremeters,Train_Ops);
    end
    %%
    if strcmp(Train_Ops.way,'LS_Opt')
        save(strcat(path_base,'【变量】\',Case_Name,'训练M矩阵结果\',KPM_Parameters_Info,'__',IDL_Peremeters.type,'.mat'),'M');
    end
end
Data_No0=Data;



%% PF Test
% %---- Generate the Testing Data
% Measure_Noise.state.input=1;% Use the noised data
% Measure_Noise.state.output=0;% Use Real data test accuracy
% fprintf('\n\nTest with the New Generated Test Data:\n');
% [Input_test,Output_test] = Generate_Tranning_OR_Testing_Data(input_num,test_input_num,num_test_input,Testing_Range,Device_Info,case_name,case_name_or_simp,Measure_Noise,TestItem);

% %-------Load in Test Data-----------
load(strcat(path_base,'【变量】\',Case_Name,'训练数据\',testing_data_file,'.mat'));
Input_Test=Data.Input;
Output_Test=Data.Output;
[Output_Test_KPMCaled] = Case_PF_Accuracy_Test(IDL_Peremeters,M_Incompelet,M_Struct,Input_Test);
fprintf('\nTest DD IDL PF， Average Error： %f   ;Maximum Error:  %f  ',mean(mean(abs(Output_Test_KPMCaled-Output_Test))),max(max(abs(Output_Test_KPMCaled-Output_Test))));



%% No0 Test
%=========================================
pos_0=find(Data.Input(:,1)==0);

Data_No0.Input(pos_0,:)=[];
IDL_Peremeters2=IDL_Peremeters;
IDL_Peremeters2.Un_RisePos=Pos_In_OutPut.Un_RisePos-length(pos_0);
IDL_Peremeters2.RisePos=Find_Wihtout_Which([1:size(Data_No0.Input,1)],Pos_In_OutPut.Un_RisePos);
IDL_Peremeters2.cent = Generate_Cent(IDL_Peremeters2,Nrbf1);
IDL_Peremeters2.type;
[M_Incompelet22,M_Struct22,M_Output22] = Cal_KMP_Total_QS(Data_No0,IDL_Peremeters2,Train_Ops);

Data_No0_Test=Data;
Data_No0_Test.Input(pos_0,:)=[];
Input_Test=Data_No0_Test.Input;
Output_Test=Data_No0_Test.Output;
[Output_Test_KPMCaled] = Case_PF_Accuracy_Test(IDL_Peremeters2,M_Incompelet22,M_Struct22,Input_Test);
fprintf('\nTest DD IDL PF， Average Error： %f   ;Maximum Error:  %f  ',mean(mean(abs(Output_Test_KPMCaled-Output_Test))),max(max(abs(Output_Test_KPMCaled-Output_Test))));
%=========================================


% [~,D1]=eig(Input_Test*Input_Test');
% min(diag(D1))
% [Lifted_Vector] = Lift_Vector_Incomplete_Total(Input_Test,IDL_Peremeters);
% Temp3=Lifted_Vector.Input_Lifted;
% % Temp3(13,:)=[];
% Temp3=Temp3(1:89,:);
% [~,D1]=eig(Temp3*Temp3');
% min(diag(D1))
% %-------Save Data-----------
% Error_PF_Test.data_Caled=Output_Test_KPMCaled;Error_PF_Test.data_standard=Output_Test;Error_PF_Test.def=abs(Output_Test_KPMCaled-Output_Test);
% path_temp=strcat(path_last,'Koopman_PF_Dispatch\存储数据\IDL等等的潮流计算误差\',Case_Name,'__',training_data_file,'\');
% mkdir(path_temp);
% save(strcat(path_temp,IDL_Peremeters.type,'.mat'),'Error_PF_Test');
% %-------Draw_PDF-----------
% Draw_PDF(mean(abs(Output_Test_KPMCaled-Output_Test)),0:0.001:0.01,1);
% Draw_PDF(max(abs(Output_Test_KPMCaled-Output_Test)),0:0.001:0.02,1);

%%-----------【Compelet Dimension Lfted (CDL) PF Test】----------------------------------
% % -------Train-----------
% cent_temp = rand(test_input_num,Nrbf1);
% cent_output = rand(size(Output,1),Nrbf2);
% [M_Compelet,M_Inverse]=Cal_M_Koopman(Data.Input,Data.Output,IDL_Peremeters.rbf_type,rbf_type_Inverse,cent_temp,cent_output);
% % -------Test-----------
% Lifted_Output_Test_Compelet = Lift_Vector_Complete(Input_Test,IDL_Peremeters.rbf_type,cent_temp);
% Output_Test_Compelet=M_Compelet*Lifted_Output_Test_Compelet;
% fprintf('\nTest DD IDL PF， Average Error： %f   ;Maximum Error:  %f  ',mean(mean(abs(Output_Test_Compelet-Output_Test))),max(max(abs(Output_Test_Compelet-Output_Test))));
% % -------Save Data-----------
% Error_PF_Test.data_Caled=Output_Test_Compelet;Error_PF_Test.data_standard=Output_Test;Error_PF_Test.def=abs(Output_Test_Compelet-Output_Test);
% path_temp=strcat(path_last,'Koopman_PF_Dispatch\存储数据\IDL等等的潮流计算误差\',Case_Name,'__',training_data_file,'\');
% mkdir(path_temp);
% % save(strcat(path_temp,'CDL','.mat'),'Error_PF_Test');
%%
Runtimes=1000;
PPeneration_Index=zeros(Runtimes,1);
PReverse_Index=zeros(Runtimes,1);
for i_time_count=1:Runtimes


        LoadinCase_Info.Case_Name=Case_Name;
        LoadinCase_Info.case_type='Massive';  %'Single'  'Massive'
        LoadinCase_Info.case_chose='case116';   %(JUST used in single test)%'case116'  'case109' % Different DG and Load Power;
        LoadinCase_Info.number=i_time_count;
        if 0%strcmp(Case_Name,'IEEE85')||strcmp(Case_Name,'IEEE141')【Position1】
            LoadinCase_Info.case_S_type='ALL_Voltage'; %'ALL_Voltage'(both below & over voltage)
        else
            LoadinCase_Info.case_S_type='OverVoltage'; %'OverVoltage'
        end
        Input_Standard=Loadin_Cases(path_last,LoadinCase_Info);  % (Generation Function  F:\科研\Koopman_PF_Dispatch\Generate_TestCases.m)
    PPeneration_Index(i_time_count,1)=sum(Input_Standard(Pos_In_OutPut.Pos_P_in_Inventer))/sum(Input_Standard(Pos_In_OutPut.Pos_P_Load));
    PReverse_Index(i_time_count,1)=sum(Input_Standard(Pos_In_OutPut.Pos_P_in_Inventer))-sum(Input_Standard(Pos_In_OutPut.Pos_P_Load));

end
fprintf(' \n  Power Peneration: max: %f   ;,min:  %f',max(PPeneration_Index),min(PPeneration_Index));
fprintf(' \n  Power Reverse: max: %f   ;,min:  %f',max(PReverse_Index),min(PReverse_Index));


%%
i_time_count=1;%Nothing;for no rolling
%% ---------------Rolling Test Accuracy-1-------------------------

Index_Weathear_Save=1;
Error_Rol.data=zeros(Runtimes,2);
Error_Rol.error_info=cell(Runtimes,1);
Error_Rol.Value_objective=zeros(Runtimes,1);
Error_Rol.OptResults.V_PQ_Cal=zeros(length(pq),Runtimes);
Error_Rol.OptResults.V_PQ_Check=zeros(length(pq),Runtimes);
Error_Rol.OptResults.Q_adj=zeros(length(Pos_In_OutPut.Pos_Q_in_Inventer),Runtimes);
Error_Rol.solve_time=zeros(Runtimes,1);
for i_time_count=start_in_1000:Runtimes

fprintf('\nSolving the No. %d snapshot',i_time_count);
%% Case in
Qno0=0;  %1 :Have data
%------Weather Generate the New Cases / Load Case ---------------------
if 0
    [Input_Standard] = Generate_Standard_Case_Vector(case_name_or_simp,Device_Info,0.7);
    LoadinCase_Info.case_chose='New Generated Case';
else
    LoadinCase_Info.Case_Name=Case_Name;
    LoadinCase_Info.case_type='Massive';  %'Single'  'Massive'
    LoadinCase_Info.case_chose='case116';   %(JUST used in single test)%'case116'  'case109' % Different DG and Load Power;
    LoadinCase_Info.number=i_time_count; 
    if 0 % strcmp(Case_Name,'IEEE85')||strcmp(Case_Name,'IEEE141')【Position2】
        LoadinCase_Info.case_S_type='ALL_Voltage'; %'ALL_Voltage'(both below & over voltage)
    else 
        LoadinCase_Info.case_S_type='OverVoltage'; %'OverVoltage' 
    end
    Input_Standard=Loadin_Cases(path_last,LoadinCase_Info);  % (Generation Function  F:\科研\Koopman_PF_Dispatch\Generate_TestCases.m)
end
%-------Weather Put the Reactive Power to 0--------------------
if Qno0
    ;
else
    Input_Standard(Pos_In_OutPut.Pos_Q_in_Inventer,:)=0;  %Weather to put 0 in Qoutput of PV
end
%--------------------------------
% A=[Input_Standard(Pos_In_OutPut.Pos_P_in_Inventer,1),Input_Standard(Pos_In_OutPut.Pos_Q_in_Inventer,1)];
% figure;bar(A);legend('光伏输出有功','光伏输出无功');
% xlabel('光伏节点序号');ylabel('功率值/MVA');title('标准算例原始值')
% set(gca,'XTicklabel',Device_Info.INclude_PV_node_str)



%% Case Basic Analyze
Input_Standard(Pos_In_OutPut.Pos_TransTab,:)=1;
Input_Standard(Pos_In_OutPut.Pos_C_Bank,:)=0;  % C bank put to 0; indeed it has no influences cos its on the referrence bus

[Output_NR,Input_NR,casename_NR,results_NR]=NR_PF_Cal(case_name,Input_Standard,Device_Info);

[Output_Test_Incompelet] = Case_PF_Accuracy_Test(IDL_Peremeters,M_Incompelet,M_Struct,Input_Standard);

% figure;plot(results_NR.bus(Pos_In_OutPut.pq,8));
% hold on;plot(Output_Test_Standard(pq-1,1));legend('NR_Method','Koopman Incomplete Lifting')
% xlabel('Bus ID'),ylabel('Voltage Magnitude(p.u.)')
% set(gcf,'unit','centimeters','position',[30,1,13,10])
% figure;plot(Input,'b')
% hold on;plot(Input_Standard,'r')
% xlabel('Inout variable'),ylabel('value');legend('The Tranning data','the standard test data')
% set(gcf,'unit','centimeters','position',[30,15,13,10])
%% Generate The Cplex Case
Constraint.PVinverter_P=Input_Standard(Pos_In_OutPut.Pos_P_in_Inventer,1)*1; %MW
Constraint.PVinverter_Q=Input_Standard(Pos_In_OutPut.Pos_Q_in_Inventer,1)*1; %MW
Constraint.P_Load=[0;Input_Standard(Pos_In_OutPut.Pos_P_Load,1)]*1;
Constraint.Q_Load=[0;Input_Standard(Pos_In_OutPut.Pos_Q_Load,1)]*1;
Constraint.INclude_PV_node=Device_Info.INclude_PV_node;
Constraint.INclude_PV_S=Device_Info.INclude_PV_S';%MVA
Constraint.Voltage.max=1.05;
Constraint.Voltage.min=0.95;
Objective_Info.Hybrid_Index=[1,0.1];
Objective_Info.Voltage_Tatget=ones(size(pq,1),1)*results_NR.gen(1,6);
Objective_Info.Opt_Variable=Opt_Variable;

switch IDL_Peremeters.type	% 1:the basic incomplete dimension lifting; 2: The incomplete dimensiong lifting with senstive variable
    case 'IDL'
        switch Opt_Variable
            case 'Reactive_Power'
                [Result_OPT_KPM,M_Matrix] = Cplex_KPM2_Opt_Only_PV_Q(Input_Standard,M_Incompelet,Pos_In_OutPut,IDL_Peremeters ,Constraint,Objective_Info);
            case 'Active_Power'
                [Result_OPT_KPM,M_Matrix] = Cplex_KPM2_Opt_Only_PV_P(Input_Standard,M_Incompelet,Pos_In_OutPut,IDL_Peremeters ,Constraint,Objective_Info);
            case 'DG_PG_Power'
                [Result_OPT_KPM,M_Matrix] = Cplex_KPM2_Opt_Only_PV_PQ(Input_Standard,M_Incompelet2,Pos_In_OutPut,cent_input_incompelet ,Constraint,Objective_Info);
        end
    case 'IDL_QuaModified'
         [Result_OPT_KPM,M_Matrix] = Cplex_KPM4_Opt_Only_PV_Q_bothPQ(Input_Standard,M_Incompelet,Pos_In_OutPut,IDL_Peremeters ,Constraint,Objective_Info);
    case 'IDLSM2'
        switch Opt_Variable
            case 'Reactive_Power'
                [Result_OPT_KPM,M_Matrix] = Cplex_KPM7_Opt_Only_PV_Q(Input_Standard,M_Incompelet,Pos_In_OutPut,IDL_Peremeters ,Constraint,Objective_Info);
        end
    case 'IDL_CubicSpace'
        switch Opt_Variable
            case 'Reactive_Power'
                [Result_OPT_KPM,M_Matrix] = Cplex_KPM8_Opt_Only_PV_Q(Input_Standard,M_Incompelet,Pos_In_OutPut,IDL_Peremeters ,Constraint,Objective_Info);
        end
end

% 【No Checking】[Result_OPT_KPM_XXX,M_Matrix_XXX] = Cplex_KPM_Opt_Total__Only_PV_Q(Input_Standard,M_Output,Pos_In_OutPut,IDL_Peremeters ,Constraint,Objective_Info);


%% ---------------Opt Results Check--------------
input_temp=Input_Standard;
% -------------KPM check-------------
input_temp([Pos_In_OutPut.Pos_P_in_Inventer,Pos_In_OutPut.Pos_Q_in_Inventer],1)=[Result_OPT_KPM.PV_P;Result_OPT_KPM.PV_Q];

[Output_Test_Standard] = Case_PF_Accuracy_Test(IDL_Peremeters,M_Incompelet,M_Struct,input_temp);
Result_OPT_KPM.Vm_PQ_KPMCheck=Output_Test_Standard(1:length(TestItem.PQ_Vm),1);
Result_OPT_KPM.PLoss_KPMCheck=Output_Test_Standard(1+length(TestItem.PQ_Vm),1);
% -------------KPM's NR PF check-------------
[~,~,~,results_NR_AfterAdj]=NR_PF_Cal(case_name,input_temp,Device_Info);
Result_OPT_KPM.Vm_PQ_NRCheck=results_NR_AfterAdj.bus(pq,8);
Result_OPT_KPM.PLoss_NRCheck=sum(abs(results_NR_AfterAdj.branch(:,14)+results_NR_AfterAdj.branch(:,16)));
% -------------KPM's NR PF SenstiveMatrix check (ONLY for DG)-------------
% % Before adjust
% [~,~,~,results_NR_BeforeAdj]=NR_PF_Cal(case_name,Input_Standard,Device_Info);
% [~,~,~,SensiM]=makeJac_20200713(results_NR_BeforeAdj.baseMVA, results_NR_BeforeAdj.bus, results_NR_BeforeAdj.branch, results_NR_BeforeAdj.gen, 1);
% SensMatrix_NR=SensiM.VM(Pos_In_OutPut.pq,Pos_In_OutPut.Bus_Num+Device_Info.INclude_PV_node);
% figure;
% subplot(2,2,1);surf(M_Matrix.M11_AfterAdj);alpha(0.25);hold on; surf(SensMatrix_NR/case_name.baseMVA);alpha(0.25);zlabel('senstive value');title('Before Adjust');
% subplot(2,2,2);surf(abs(M_Matrix.M11_AfterAdj-SensMatrix_NR/case_name.baseMVA));zlabel('error');title('Before Adjust');
% fprintf('\n 【After Adjust】The error of senstiven matrix(only the DG bus),Average Error: %f ;Maximum Error: %f ;',mean(mean(abs(M_Matrix.M11_AfterAdj-SensMatrix_NR/case_name.baseMVA))),max(max(abs(M_Matrix.M11_AfterAdj-SensMatrix_NR/case_name.baseMVA))));
% 
% % After adjust
% [~,~,~,SensiM]=makeJac_20200713(results_NR_AfterAdj.baseMVA, results_NR_AfterAdj.bus, results_NR_AfterAdj.branch, results_NR_AfterAdj.gen, 1);
% SensMatrix_NR=SensiM.VM(Pos_In_OutPut.pq,Pos_In_OutPut.Bus_Num+Device_Info.INclude_PV_node);
% % figure;
% subplot(2,2,3);surf(M_Matrix.M11);alpha(0.25);hold on; surf(SensMatrix_NR/case_name.baseMVA);alpha(0.25);zlabel('senstive value');title('After Adjust');
% subplot(2,2,4);surf(abs(M_Matrix.M11-SensMatrix_NR/case_name.baseMVA));zlabel('error');title('After Adjust');
% fprintf('\n 【Before Adjust】The error of senstiven matrix(only the DG bus),Average Error: %f ;Maximum Error: %f ;',mean(mean(abs(M_Matrix.M11-SensMatrix_NR/case_name.baseMVA))),max(max(abs(M_Matrix.M11-SensMatrix_NR/case_name.baseMVA))));
%% ---------------Rolling Test Accuracy--2------------------------
fprintf('\n Adjustment Error： Average: %f ;  Maximum: %f',mean(abs([results_NR_AfterAdj.bus(1,8);Result_OPT_KPM.Vm_PQ]-results_NR_AfterAdj.bus(:,8))),...
    max(abs([results_NR_AfterAdj.bus(1,8);Result_OPT_KPM.Vm_PQ]-results_NR_AfterAdj.bus(:,8))));
Error_Rol.data(i_time_count,1)=mean(abs([results_NR_AfterAdj.bus(1,8);Result_OPT_KPM.Vm_PQ]-results_NR_AfterAdj.bus(:,8)));
Error_Rol.data(i_time_count,2)=max(abs([results_NR_AfterAdj.bus(1,8);Result_OPT_KPM.Vm_PQ]-results_NR_AfterAdj.bus(:,8)));
Error_Rol.error_info{i_time_count,1}=Result_OPT_KPM.Opt_Info;
Error_Rol.Value_objective(i_time_count,1)=Result_OPT_KPM.objective;
Error_Rol.OptResults.V_PQ_Cal(:,i_time_count)=Result_OPT_KPM.Vm_PQ;
Error_Rol.OptResults.V_PQ_Check(:,i_time_count)=Result_OPT_KPM.Vm_PQ_NRCheck;
Error_Rol.OptResults.Q_adj(:,i_time_count)=Result_OPT_KPM.PV_Q_Adj;
Error_Rol.OptResults.P_adj(:,i_time_count)=Result_OPT_KPM.PV_P_Adj;
% Error_Rol.VDR_Check=sum(abs(Objective_Info.Voltage_Tatget-Result_OPT_KPM.Vm_PQ_NRCheck))/length(pq);
Error_Rol.solve_time(i_time_count,1)=Result_OPT_KPM.solvertime;
    %% Single Save

if Index_Weathear_Save
    Error_Rol.ID_count=i_time_count;
    path_temp=strcat(path_last,'Koopman_PF_Dispatch\存储数据\大规模算例计算误差积累分布函数计算结果\',KPM_Parameters_Info,'\Temp_Single\');
    mkdir(path_temp)
    save(strcat(path_last,'Koopman_PF_Dispatch\存储数据\大规模算例计算误差积累分布函数计算结果\',KPM_Parameters_Info,'\Temp_Single\',case_obj_info,'.mat'),'Error_Rol')
end
end

Index_Solveable_Index_temp = Judge_Solveable_Cell(Error_Rol.error_info);
fprintf('\n AE:%f ,ME: %f ',mean(Error_Rol.data(Index_Solveable_Index_temp,1)),max(Error_Rol.data(Index_Solveable_Index_temp,2)))
        
temp=sort(Error_Rol.data(:,1));Error_Rol.MeanError_p1=temp(ceil(length(temp)*0.01));Error_Rol.MeanError_p50=temp(ceil(length(temp)*0.5));Error_Rol.MeanError_p99=temp(ceil(length(temp)*0.99));
temp=sort(Error_Rol.data(:,2));Error_Rol.MaxError_p1=temp(ceil(length(temp)*0.01));Error_Rol.MaxError_p50=temp(ceil(length(temp)*0.5));Error_Rol.MaxError_p99=temp(ceil(length(temp)*0.99));
figure;subplot(1,2,1);cdfplot(Error_Rol.data(:,1));title('Average Error PDF');subplot(1,2,2);cdfplot(Error_Rol.data(:,2));title('Maximum Error PDF');
if Index_Weathear_Save
    path_temp=strcat(path_last,'Koopman_PF_Dispatch\存储数据\大规模算例计算误差积累分布函数计算结果\',KPM_Parameters_Info,'\');
    mkdir(path_temp)
    save(strcat(path_last,'Koopman_PF_Dispatch\存储数据\大规模算例计算误差积累分布函数计算结果\',KPM_Parameters_Info,'\',case_obj_info,'.mat'),'Error_Rol')
end
% END %---------------Rolling Test Accuracy--------------------------


 
% %% Drawn Results Figure & Case Comparison
% No_Check_and_Comparison=0;
% if No_Check_and_Comparison
%     fprintf('\n不做算例对比\n');
%     figure;plot(results_NR.bus(Pos_In_OutPut.pq,8));hold on;plot(Result_OPT_KPM.Vm_PQ);hold on;plot( Result_OPT_KPM.Vm_PQ_NRCheck);
%     title('case116');xlabel('节点编号');ylabel('节点电压p.u.');legend('优化前','优化后','优化后校验');set(gcf,'unit','centimeters','position',[1,1,13,10])
%     switch Opt_Variable
%         %  Pos_In_OutPut.Pos_Q_in_Inventer  Pos_In_OutPut.Pos_P_in_Inventer  [Pos_In_OutPut.Pos_P_in_Inventer;Pos_In_OutPut.Pos_Q_in_Inventer]
%         case 'Reactive_Power'
%             figure;bar(Result_OPT_KPM.PV_Q);xlabel('Bus ID');ylabel('PV Reactive power injection /MVar');title('Results Comparison ');set(gcf,'unit','centimeters','position',[15,1,13,10])
%             figure;bar(Result_OPT_KPM.PV_Q_Adj);xlabel('PV inventer ID');ylabel('Reactive Power / MVar');title('Adjustment （Reactive Power of Photovoltaic）');set(gcf,'unit','centimeters','position',[15,15,13,10])
%         case 'Active_Power'
%             figure;bar(Result_OPT_KPM.PV_P);xlabel('Bus ID');ylabel('PV Active power injection /MVar');title('Results Comparison ');set(gcf,'unit','centimeters','position',[15,1,13,10])
%             figure;bar(Result_OPT_KPM.PV_P_Adj);xlabel('PV inventer ID');ylabel('Active Power / MVar');title('Adjustment （Active Power of Photovoltaic）');set(gcf,'unit','centimeters','position',[15,15,13,10])
%     end
% else
%     fprintf('\n做算例对比\n');
%     %% ------------DLPF Result and NR Check-----------------
%     switch Case_Name
%     case 'IEEE33'
%         path_storge=strcat(path_last,'Koopman_PF_Dispatch\存储数据\主函数输出数据\Case33\');
%     case 'IEEE69'
%         path_storge=strcat(path_last,'Koopman_PF_Dispatch\存储数据\主函数输出数据\Case69\');
%         
%     end
%     [temp_load_data,temp_load_data_DDLPF] = Loadin_Data(path_storge,Objective_Info.type,LoadinCase_Info.case_chose,Qno0);%temp_load_data:Model Based Method  temp_load_data_DDLPF:Data Based Method
%     %% Data Transfer 1 Model Based Method 
%     [CPLEX_MBM_Results] = Loaded_Data_Preprocess(temp_load_data,case_name,Input_Standard,Pos_In_OutPut,Device_Info);
%      %% Data Transfer 2 Data Based Method
%     [CPLEX_DDM_Results] = Loaded_Data_Preprocess(temp_load_data_DDLPF,case_name,Input_Standard,Pos_In_OutPut,Device_Info);
%     %% Results Comparison
%     Results_Comparison_KPM_Opt(results_NR,Output_Test_Standard,Result_OPT_KPM,CPLEX_MBM_Results,CPLEX_DDM_Results,Objective_Info,Pos_In_OutPut,LoadinCase_Info.case_chose);    
%     %%
%     switch Opt_Variable
%         %  Pos_In_OutPut.Pos_Q_in_Inventer  Pos_In_OutPut.Pos_P_in_Inventer  [Pos_In_OutPut.Pos_P_in_Inventer;Pos_In_OutPut.Pos_Q_in_Inventer]
%         case 'Reactive_Power'
%             fprintf('\nKMP优化后调节量：%f , DLPF/SOCP 优化后调节量:%f , DDLPF优化后调节量:%f',sum(abs(Result_OPT_KPM.PV_Q_Adj)),sum(abs(CPLEX_MBM_Results.PV_Q_Adj)),sum(abs(CPLEX_DDM_Results.PV_Q_Adj)));
%         case 'Active_Power'
%             fprintf('\nKMP优化后调节量：%f , DLPF/SOCP 优化后调节量:%f , DDLPF优化后调节量:%f',sum(abs(Result_OPT_KPM.PV_P_Adj)),sum(abs(CPLEX_MBM_Results.PV_P_Adj)),sum(abs(CPLEX_DDM_Results.PV_Q_Adj)));
%     end
% end
% 








