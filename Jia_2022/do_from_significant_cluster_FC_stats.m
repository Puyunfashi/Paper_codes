% This script does stats based on FC analysis from one cluster 
% Xiao Chen 220429
% chenxiaophd@gmail.com

%% initialization
clear;clc;

data_dir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Preprocessing';
work_dir = ['/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/', ...
                                    'Regional_Metrics_HeadMotion_Updated/corr_analysis_V4'];
TemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/CX_software/DPABI-master/DPABISurf/SurfTemplates';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end

%% select participants 
demographic_info = readtable(['/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/', ...
                                                                                        'Demographic_Info_V3.txt'],'Delimiter','\t');

sub_list = demographic_info.Serial_Number;
Dx =  demographic_info.Dx;
Age = demographic_info.Age;
Sex = demographic_info.Sex;
Sex = Sex -1;
Edu = demographic_info.Edu;

% some subjects do not have rum and dis data
% read in subjects who have rum/dis data
sub_list_rum_dis = importdata(['/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/', ...
                                                    'Preprocessing/Task/sub_list_task_merged20220118.txt']);
% exclude those who don't have rum and dis data
WantedSubMatrix = ismember(sub_list, sub_list_rum_dis);
WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);

% rule out subjects with missing data, xx = 999
WantedSubMatrix = ones(length(sub_list),1);
WantedSubMatrix(Dx == 999) = 0;
WantedSubMatrix(Age == 999) = 0;
WantedSubMatrix(Sex == 999) = 0;
WantedSubMatrix(Edu == 999) = 0;

WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);

% Extract head motion
% 1: Rum, 2: Dis
HeadMotion = [];
for i=1:length(sub_list)
    Temp=load([data_dir,'/Task/RealignParameter/sub-',sub_list{i},'/FD_Jenkinson_sub-',sub_list{i},'.txt']);
    HeadMotion(i,1)=mean(Temp);
     Temp=load([data_dir,'/Task/RealignParameter/sub-',sub_list{i},'/S2_FD_Jenkinson_sub-',sub_list{i},'.txt']);
    HeadMotion(i,2)=mean(Temp);
end
% rule out those headmotion > 0.2
WantedSubMatrix = ones(length(sub_list),1);
WantedSubMatrix(HeadMotion(:,1)>0.2) = 0;
WantedSubMatrix(HeadMotion(:,2)>0.2) = 0;
WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);
HeadMotion = HeadMotion(WantedSubIndex,:);

% Sub026 and Sub051's phase coding directions are different
% Sub005 is left handed
WantedSubMatrix = ones(length(sub_list),1);
WantedSubMatrix(strcmp(sub_list,  'Sub026')) = 0;
WantedSubMatrix(strcmp(sub_list,  'Sub051')) = 0;
WantedSubMatrix(strcmp(sub_list,  'Sub005')) = 0;
WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);
HeadMotion = HeadMotion(WantedSubIndex,:);

%% Arrange Files 
% Arrange FC Files
TargetDirRum = [work_dir,'/FC_Maps/MDD_Rum'];
if ~exist(TargetDirRum); mkdir(TargetDirRum); end
mkdir([TargetDirRum,'/lh']); mkdir([TargetDirRum,'/rh']); mkdir([TargetDirRum,'/subcorticol']); 
TargetDirDis = [work_dir,'/FC_Maps/MDD_Dis'];
mkdir([TargetDirDis,'/lh']); mkdir([TargetDirDis,'/rh']); mkdir([TargetDirDis,'/subcorticol']); 
if ~exist(TargetDirDis); mkdir(TargetDirDis); end

for i = 1:length(find(Dx == 1))
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/Results/FunSurfLH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirRum,'/lh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/Results/FunSurfRH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirRum,'/rh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/zFC_sub-',sub_list{i},'.nii'];
    TargetFile = [TargetDirRum,'/subcorticol/zFC_sub-',sub_list{i},'.nii'];
    copyfile(OriginFile,TargetFile);
    
     OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/S2_Results/FunSurfLH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirDis,'/lh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/S2_Results/FunSurfRH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirDis,'/rh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/S2_Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/zFC_sub-',sub_list{i},'.nii'];
    TargetFile = [TargetDirDis,'/subcorticol/zFC_sub-',sub_list{i},'.nii'];
    copyfile(OriginFile,TargetFile);
end

TargetDirRum = [work_dir,'/FC_Maps/HC_Rum'];
if ~exist(TargetDirRum); mkdir(TargetDirRum); end
mkdir([TargetDirRum,'/lh']); mkdir([TargetDirRum,'/rh']); mkdir([TargetDirRum,'/subcorticol']); 
TargetDirDis = [work_dir,'/FC_Maps/HC_Dis'];
mkdir([TargetDirDis,'/lh']); mkdir([TargetDirDis,'/rh']); mkdir([TargetDirDis,'/subcorticol']); 
if ~exist(TargetDirDis); mkdir(TargetDirDis); end

for i = length(find(Dx == 1))+1:length(sub_list)
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/Results/FunSurfLH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirRum,'/lh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/Results/FunSurfRH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirRum,'/rh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/zFC_sub-',sub_list{i},'.nii'];
    TargetFile = [TargetDirRum,'/subcorticol/zFC_sub-',sub_list{i},'.nii'];
    copyfile(OriginFile,TargetFile);
    
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/S2_Results/FunSurfLH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirDis,'/lh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/S2_Results/FunSurfRH/FC_SeedSurfLHSurfRHVolu_FunSurfWCFS/zFC_sub-',sub_list{i},'.func.gii'];
    TargetFile = [TargetDirDis,'/rh/zFC_sub-',sub_list{i},'.func.gii'];
    copyfile(OriginFile,TargetFile);
    OriginFile = []; TargetFile = [];
    OriginFile = [data_dir,'/Task/S2_Results/FunVolu/FC_SeedSurfLHSurfRHVolu_FunVoluWCFS/zFC_sub-',sub_list{i},'.nii'];
    TargetFile = [TargetDirDis,'/subcorticol/zFC_sub-',sub_list{i},'.nii'];
    copyfile(OriginFile,TargetFile);
end

%% Arrange DC Maps
HemisphereSet = {'LH','RH'};MeasureSet = {'DegreeCentrality'};
SuffixSet = {'_Bilateral_PositiveWeightedSumBrain'};PipesuffixSet = {'_FunSurfWCF'};
HemisphereNameSet = {'lh','rh'};

%Rumination
%MDD
for i = 1:length(find(Dx == 1))
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                                MeasureSet{iMeasure},PipesuffixSet{iMeasure},'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                                '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/DC_Maps/MDD_Rum/',HemisphereNameSet{iHem}];
            if ~exist(TargetDir, 'dir'); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%HC
for i = length(find(Dx == 1))+1:length(sub_list)
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                                MeasureSet{iMeasure},PipesuffixSet{iMeasure},'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                                '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/DC_Maps/HC_Rum/',HemisphereNameSet{iHem}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%Distraction
%MDD
for i = 1:length(find(Dx == 1))
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/S2_ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                             MeasureSet{iMeasure},PipesuffixSet{iMeasure},'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                             '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/DC_Maps/MDD_Dis/',HemisphereNameSet{iHem}];
            if ~exist(TargetDir, 'dir'); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%HC
for i = length(find(Dx == 1))+1:length(sub_list)
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/S2_ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                                MeasureSet{iMeasure},PipesuffixSet{iMeasure},'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                                '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/DC_Maps/HC_Dis/',HemisphereNameSet{iHem}];
            if ~exist(TargetDir, 'dir'); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%% do mixed effect analysis on surface
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_lh_white.surf.gii'];
PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_lh_white_avg.area.gii'];
MaskFile =  [TemplateDir,'/fsaverage5_lh_cortex.label.gii'];

DependentDir{1,1} =  [work_dir,'/MDD_Rum/lh'];
DependentDir{2,1} = [work_dir,'/MDD_Dis/lh'];
DependentDir{3,1} = [work_dir,'/HC_Rum/lh'];
DependentDir{4,1} = [work_dir,'/HC_Dis/lh'];

OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

OutputDir = [work_dir,'/stats/lh'];
if ~exist(OutputDir); mkdir(OutputDir); end
OutputName = [OutputDir,'/FC_stats'];
y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);

PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_rh_white.surf.gii'];
PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_rh_white_avg.area.gii'];
MaskFile =  [TemplateDir,'/fsaverage5_rh_cortex.label.gii'];

DependentDir{1,1} =  [work_dir,'/MDD_Rum/rh'];
DependentDir{2,1} = [work_dir,'/MDD_Dis/rh'];
DependentDir{3,1} = [work_dir,'/HC_Rum/rh'];
DependentDir{4,1} = [work_dir,'/HC_Dis/rh'];

OutputDir = [work_dir,'/FC_analysis/stats/rh'];
if ~exist(OutputDir); mkdir(OutputDir); end
OutputName = [OutputDir,'/FC_stats'];
y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);


PALMSettings = [];
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

MaskFile =  [data_dir,'/Task/Masks/AllResampled_BrainMask_05_91x109x91.nii'];

DependentDir{1,1} =  [work_dir,'/MDD_Rum/subcorticol'];
DependentDir{2,1} = [work_dir,'/MDD_Dis/subcorticol'];
DependentDir{3,1} = [work_dir,'/HC_Rum/subcorticol'];
DependentDir{4,1} = [work_dir,'/HC_Dis/subcorticol'];

OutputDir = [work_dir,'/FC_analysis/stats/subcorticol'];
if ~exist(OutputDir); mkdir(OutputDir); end
OutputName = [OutputDir,'/FC_stats'];
y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);

%%  do 7 network analysis with Schaefer 400 parcelations
TemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/Matrix_Extraction';
dpabiSurfTemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/CX_software/DPABI-master/DPABISurf/SurfTemplates';
HemisphereNameSet = {'lh', 'rh'};
ConditionSet = {'MDD_Rum', 'MDD_Dis', 'HC_Rum', 'HC_Dis'};

%% extract network level FC values
TopOutputDir = [work_dir, '/FC_network_extracted'];
if ~exist(TopOutputDir); mkdir(TopOutputDir); end
% surface
IsMultipleLabel = 1;
GHeader = [];
CUTNUMBER = 1;
ROIDef = {};
for iCondition = 1:length(ConditionSet)
    for iHem = 1:2
        ROIDef{1} = [TemplateDir, '/fsaverage5_', HemisphereNameSet{iHem}, ...
                                                            '_Schaefer2018_400Parcels_17Networks_order.label.gii'];
        AllVolume = [work_dir, '/FC_Maps/', ConditionSet{iCondition}, '/', HemisphereNameSet{iHem}];
        OutputName = [TopOutputDir, '/', ConditionSet{iCondition}, '_', HemisphereNameSet{iHem}];
        AMaskFilename = [dpabiSurfTemplateDir, '/fsaverage5_', HemisphereNameSet{iHem}, '_cortex.label.gii'];
        [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, AMaskFilename, ...
                                                                    IsMultipleLabel, GHeader, CUTNUMBER);
    end
end

% subcortical
IsNeedDetrend = 0;
TemporalMask = [];
for iCondition = 1:length(ConditionSet)
    ROIDef{1} = [TemplateDir, '/Tian_Subcortex_S4_3T_2009cAsym.nii'];
    AllVolume = [work_dir, '/FC_Maps/', ConditionSet{iCondition}, '/subcorticol'];
    OutputName = [TopOutputDir, '/', ConditionSet{iCondition}, '_subcortical'];
    MaskData = [data_dir,'/Task/Masks/AllResampled_BrainMask_05_91x109x91.nii'];
    [ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, MaskData, ...
                                            IsMultipleLabel, IsNeedDetrend, [], [], TemporalMask, [], [], [], CUTNUMBER);
end

%% get network FC
load(['/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/Matrix_Extraction/', ...
                'DPABISurf_Schaefer2018_400_17Networks_Tian2020_54_Info438_V3.mat']);

% network_FC: 1 x 4 cells, 4 conditions in ConditionSet,
% each cell has n x 10 matrixs, 10 networks defined as the abovementrioned
% .mat file
TopOutputDir = [work_dir, '/FC_network_extracted'];
network_FC = {};
for iCondition = 1:length(ConditionSet)
    load([TopOutputDir,'/ROISignals_', ConditionSet{iCondition}, '_lh.mat']);
    ROISignals_lh = ROISignals;
    load([TopOutputDir,'/ROISignals_', ConditionSet{iCondition}, '_rh.mat']);
    ROISignals_rh = ROISignals;
    load([TopOutputDir,'/ROISignals_', ConditionSet{iCondition}, '_subcortical.mat']);
    ROISignals_subcortical = ROISignals;
    ROISignals_full{iCondition} = [ROISignals_lh,ROISignals_rh,ROISignals_subcortical];
    ROISignals_full_438{iCondition} = ROISignals_full{iCondition}(:, ROIIndex_181_Schafer_Tian);
    network_label = DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork';
    for iNetwork = unique(network_label)
        network_FC{iCondition}(:,iNetwork) = mean(ROISignals_full_438{iCondition}(:, find(network_label == iNetwork)), 2);
    end
end

%% export .txt files for plotting with ggplot2
for iNetwork = 1:size(network_FC{1},2)
    txt_matrix = [];
    for iCondition = 1:4
        temp = network_FC{iCondition}(:, iNetwork);
        temp = [temp,ones(length(temp),1).*iCondition];
        txt_matrix = [txt_matrix; temp];
    end
    save([work_dir,'/',YeoSCNetwork_Label{iNetwork,2},'.txt'], 'txt_matrix', '-ASCII', '-tabs');
end

%% do mixed effect analysis
MDD_Rum_n = size(network_FC{1},1);
MDD_Dis_n = size(network_FC{2},1);
HC_Rum_n = size(network_FC{3},1);
HC_Dis_n = size(network_FC{4},1);

BetweenSubjectFactor1 = ones(MDD_Rum_n+MDD_Dis_n,1);
BetweenSubjectFactor2 = -1*ones(HC_Rum_n+HC_Dis_n,1);
BetweenSubjectFactor = [BetweenSubjectFactor1;BetweenSubjectFactor2];
WithinSubjectFactor = [ones(MDD_Rum_n,1);-1*ones(MDD_Dis_n,1);ones(HC_Rum_n,1);-1*ones(HC_Dis_n,1)];
nSubject = [MDD_Rum_n MDD_Dis_n HC_Rum_n HC_Dis_n];

OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

SubjectRegressorsAll = [];
for i=1:4
    SubjectRegressorsAll=[SubjectRegressorsAll;(ceil(i/2)-1)*10000+([1:nSubject(i)]')];
end
SubIndex=unique(SubjectRegressorsAll);
nSub=length(SubIndex);
SubjectRegressors=[];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(size(SubjectRegressorsAll));
    SubjectRegressors(SubjectRegressorsAll==SubIndex(i),i) = 1;
end
Interaction = WithinSubjectFactor.*BetweenSubjectFactor;
HeadMotion_Regressor = [OtherCovariates{1};OtherCovariates{2};OtherCovariates{3};OtherCovariates{4}];
AllCov = [WithinSubjectFactor,Interaction,SubjectRegressors,HeadMotion_Regressor];
Contrast = zeros(1,size(AllCov,2));
TF_Flag = 'T';

% condition effect
Contrast(1)=1;
for iNetwork = unique(network_label)
    y = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,TF_Flag);
    Stats_Mixed_Condition{iNetwork,1} = YeoSCNetwork_Label{iNetwork,2};
    Stats_Mixed_Condition{iNetwork,2} = TF_ForContrast;
    Stats_Mixed_Condition{iNetwork,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
    Stats_Mixed_Condition{iNetwork,4} = 2*tcdf(-abs(TF_ForContrast),nSub-2)*9;
    Stats_Mixed_Condition{iNetwork,5} = Cohen_f2;
end

% interaction effect
Contrast(1)=0;
Contrast(2)=1;
for iNetwork = unique(network_label)
    y = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,TF_Flag);
    Stats_Mixed_Interaction{iNetwork,1} = YeoSCNetwork_Label{iNetwork,2};
    Stats_Mixed_Interaction{iNetwork,2} = TF_ForContrast;
    Stats_Mixed_Interaction{iNetwork,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
    Stats_Mixed_Interaction{iNetwork,4} = 2*tcdf(-abs(TF_ForContrast),nSub-2)*9;
    Stats_Mixed_Interaction{iNetwork,5} = Cohen_f2;
    p_raw(iNetwork) =  2*tcdf(-abs(TF_ForContrast),nSub-2);
end
% fdr correction
p_raw = p_raw(1:9);
p_fdr = mafdr(p_raw, 'BHFDR',true);


%% do post-hoc comparisons
% a structure stats_posthoc_networkFC with 4 fields: 
% group_effect_rum; group_effect_dis;
% condition_effect_MDD; condition_effect_HC

for iNetwork = 1:size(network_FC{1},2)
    [~,p,~,stats] = ttest(network_FC{1}(:,iNetwork), network_FC{2}(:,iNetwork));
    stats_posthoc_networkFC.condition_effect_MDD(iNetwork,1) = stats.tstat;
    stats_posthoc_networkFC.condition_effect_MDD(iNetwork,2) = p;
    
    [~,p,~,stats] = ttest(network_FC{3}(:,iNetwork), network_FC{4}(:,iNetwork));
    stats_posthoc_networkFC.condition_effect_HC(iNetwork,1) = stats.tstat;
    stats_posthoc_networkFC.condition_effect_HC(iNetwork,2) = p;
    
    [~,p,~,stats] = ttest2(network_FC{1}(:,iNetwork), network_FC{3}(:,iNetwork));
    stats_posthoc_networkFC.group_effect_rum(iNetwork,1) = stats.tstat;
    stats_posthoc_networkFC.group_effect_rum(iNetwork,2) = p;
    
    [~,p,~,stats] = ttest2(network_FC{2}(:,iNetwork), network_FC{4}(:,iNetwork));
    stats_posthoc_networkFC.group_effect_dis(iNetwork,1) = stats.tstat;
    stats_posthoc_networkFC.group_effect_dis(iNetwork,2) = p;
end

%% plot
group_label = {};
for i = 1:sum(nSubject)
    if i<=nSubject(1)
            group_label{i,1} = 'MDD Rum';
    elseif (i > nSubject(1)) && (i <= nSubject(1) + nSubject(2))
            group_label{i,1} = 'MDD Dis';
    elseif (i > nSubject(1) + nSubject(2)) && (i <= nSubject(1) + nSubject(2) + nSubject(3))
            group_label{i,1} = 'HC Rum';
    else
            group_label{i,1} = 'HC Dis';
    end
end

addpath(genpath('/mnt/Data2/RfMRILab/Chenxiao/CX_software/Violinplot-Matlab-master'));

figure('DefaultTextFontName','Helvetica','DefaultAxesFontName','Helvetica');
set(gcf,'position',[150, 150, 800, 7000]);

group_order  = {'MDD Rum', 'MDD Dis', 'HC Rum', 'HC Dis'};
for iNetwork = unique(network_label)
    subplot(5,2,iNetwork)
    data = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
    violinplot(data, group_label, 'GroupOrder', group_order);
    if Stats_Mixed_Interaction{iNetwork,4} < 0.05
        title([YeoSCNetwork_Label{iNetwork,2}, '*']);
    else
        title(YeoSCNetwork_Label{iNetwork,2});
    end
end

print(gcf,[work_dir,'/FC_fromDCCluster_network.jpg'],'-djpeg','-r600');
close all;

