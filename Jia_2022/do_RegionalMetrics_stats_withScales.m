% Some subjects' scale data are missing, so all the analyses regarding
% scales were listed here, with some MDD patients excluded.
% Written by Xiao Chen 220922
% chenxiaophd@gmail.com

%% initialization
clear;clc;

data_dir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Preprocessing';
work_dir = ['/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/', ...
                                    'Regional_Metrics_HeadMotion_Updated/scale_corr_analysis_V8'];
TemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/CX_software/DPABI-master/DPABISurf/SurfTemplates';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end

%% select participants 
demographic_info = readtable(['/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/', ...
                                                                                        'Demographic_Info_V5.txt'],'Delimiter','\t');

sub_list = demographic_info.Serial_Number;
Dx =  demographic_info.Dx;
Age = demographic_info.Age;
Sex = demographic_info.Sex;
Sex = Sex -1;
Edu = demographic_info.Edu;

% % to select all scales and blood metrics
% scale_score = [demographic_info.Rumination, demographic_info.Reflection, ...
%                         demographic_info.Brooding, ...
%                         demographic_info.SSI_C, demographic_info.SSI_W, ...
%                         demographic_info.HAMD_Sum, demographic_info.HAMA_Sum, ...
%                         demographic_info.Insomnia, demographic_info.Emotional_Abuse, ...
%                         demographic_info.Physical_Abuse, demographic_info.Sexual_Abuse, ...
%                         demographic_info.Emotional_Neglect, demographic_info.Physical_Neglect, ...
%                         demographic_info.Childhood_Trauma, demographic_info.Day_Dream, ...
%                         demographic_info.Neuroticism, demographic_info.Conscientiousness, ...
%                         demographic_info.Openness, demographic_info.Agreeableness, ...
%                         demographic_info.Extraversion, demographic_info.BDI, ...
%                         demographic_info.Negative_Events, -demographic_info.Positive_Events, ...
%                         demographic_info.All_Events, ...
%                         demographic_info.OGM_N, demographic_info.OGM_P, ...
%                         demographic_info.Family_Dynamic, demographic_info.Family_Atmosphere, ...
%                         demographic_info.Individulization, demographic_info.System_Logic, ...
%                         demographic_info.Diease_View, demographic_info.Family_Function, ...
%                         demographic_info.Problem_Solve, demographic_info.Communication, ...
%                         demographic_info.Charactor, demographic_info.Emotional_Reaction, ...
%                         demographic_info.Emontional_Interfere, demographic_info.Behavioral_Control, ...
%                         demographic_info.Function_Full, ...
%                         demographic_info.IL_1beta, demographic_info.IL_2, demographic_info.IL_6, ...
%                         demographic_info.IL_8, demographic_info.IL_18, demographic_info.TNFalpha, ...
%                         demographic_info.IFNy, demographic_info.GDNF, demographic_info.CRP
%                         ];
%                     
% scale_name_set = {  'Rumination', 'Reflection', 'Brooding',  ...
%                             'SSI_C', 'SSI_W', ...
%                             'HAMD', 'HAMA', ...
%                             'Insomnia', 'Emotional_Abuse', 'Physical_Abuse', ...
%                             'Sexual_Abuse', 'Emotional_Neglect', 'Physical_Neglect', ...
%                             'Childhood_Trauma', 'Day_Dream', 'Neuroticism', ...
%                             'Conscientiousness', 'Openness', 'Agreeableness', 'Extraversion', 'BDI' ...
%                             'Negative_Events', 'Positive_Events', 'All_Events', ...
%                             'OGM_N', 'OGM_P', 'Family_Dynamic', 'Family_Atmosphere', ...
%                             'Individulization', 'System_Logic', 'Diease_View', 'Family_Function', ...
%                             'Problem_Solve', 'Communication', 'Charactor', 'Emotional_Reaction', ...
%                             'Emontional_Interfere', 'Behavioral_Control', 'Function_Full', ...
%                             'IL_1beta','IL_2','IL_6','IL_8','IL_18','TNFalpha', 'IFNy', 'GDNF', 'CRP'
%                            };
                       
% to select some interested scales: rumination, OGM
scale_score = [demographic_info.Rumination, demographic_info.Reflection, ...
                      demographic_info.Brooding, ...
                      demographic_info.HAMD_Sum, ...
                      demographic_info.BDI, ...
                      demographic_info.OGM_N, demographic_info.OGM_P, ...
                      demographic_info.OGM_N + demographic_info.OGM_P, ...
                       ];
                   
scale_name_set = { 'Rumination', 'Reflection', 'Brooding', ...
                              'HAMD', 'BDI', ...
                              'OGM_N', 'OGM_P', 'OGM_full'
                              };

% % to select some interested scales: family, OGM
% scale_score = [demographic_info.OGM_N, demographic_info.OGM_P, ...
%                       demographic_info.OGM_N + demographic_info.OGM_P, ...
%                       demographic_info.Family_Dynamic, demographic_info.Family_Atmosphere, ...
%                      demographic_info.Individulization, demographic_info.System_Logic, ...
%                      demographic_info.Diease_View, demographic_info.Family_Function, ...
%                      demographic_info.Problem_Solve, demographic_info.Communication, ...
%                      demographic_info.Charactor, demographic_info.Emotional_Reaction, ...
%                      demographic_info.Emontional_Interfere, demographic_info.Behavioral_Control, ...
%                      demographic_info.Function_Full
%                        ];
%                    
% scale_name_set = { 'OGM_N', 'OGM_P', 'OGM_full', ...
%                               'Family_Dynamic', 'Family_Atmosphere', ...
%                             'Individulization', 'System_Logic', 'Diease_View', 'Family_Function', ...
%                             'Problem_Solve', 'Communication', 'Charactor', 'Emotional_Reaction', ...
%                             'Emontional_Interfere', 'Behavioral_Control', 'Function_Full'
%                                };

% % include only rumination scales
% scale_score = [demographic_info.Rumination, demographic_info.Reflection, ...
%                       demographic_info.Brooding
%                        ];
%                    
% scale_name_set = { 'Rumination', 'Reflection', 'Brooding'
%                               };

                         
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

scale_score = scale_score(WantedSubIndex,:);

% rule out subjects with missing data, xx = 999
WantedSubMatrix = ones(length(sub_list),1);
WantedSubMatrix(Dx == 999) = 0;
WantedSubMatrix(Age == 999) = 0;
WantedSubMatrix(Sex == 999) = 0;
WantedSubMatrix(Edu == 999) = 0;

% for scale score
[row, col] = ind2sub(size(scale_score),find(scale_score == 999));
WantedSubMatrix(unique(row)) = 0;

% rule out subjects who don't need to fill out the suicidal scale, x = 888
scale_score(scale_score == 888) = 0;

WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);

scale_score = scale_score(WantedSubIndex,:);

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

scale_score = scale_score(WantedSubIndex,:);

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

scale_score = scale_score(WantedSubIndex,:);

%% creat a heatmap showing the correlation of scale scores
R = corrcoef(scale_score(Dx==1, :));
heatmap(scale_name_set, scale_name_set, R, 'colormap', jet);

R = corrcoef(scale_score(Dx==2, :));
heatmap(scale_name_set, scale_name_set, R, 'colormap', jet);

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

%%  extract Cluster's DC values
TopOutputDir = [work_dir, '/DC_Interaction_Cluster'];
if ~exist(TopOutputDir, 'dir'); mkdir(TopOutputDir); end
% surface
IsMultipleLabel = 0;
GHeader = [];
CUTNUMBER = 1;
ROIDef = {};
for iCondition = 1:length(ConditionSet)
    for iHem = 1:2
        ROIDef{1} = [work_dir, '/DC_interaction_Cluster_Mask.gii'];
        AllVolume = [work_dir, '/DC_Maps/', ConditionSet{iCondition}, '/', HemisphereNameSet{iHem}];
        OutputName = [TopOutputDir, '/', ConditionSet{iCondition}, '_', HemisphereNameSet{iHem}];
        AMaskFilename = [dpabiSurfTemplateDir, '/fsaverage5_', HemisphereNameSet{iHem}, '_cortex.label.gii'];
        [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, AMaskFilename, ...
                                                                    IsMultipleLabel, GHeader, CUTNUMBER);
    end
end

%% construct MRI Phenotype
% a 1 x 8 cell, each cell contains a n x 11 matrix: column 1 - 10, 10
% networks in YeoSCNetwork_Label, 11 is the DC interaction cluster value
% cell 1-4: MDD_Rum, MDD_Dis, HC_Rum, HC_Dis
% cell 5-6: MDD_Rum-MDD_Dis, HC_Rum-HC_Dis
% cell 7-8: (MDD_Rum-MDD_Dis)/MDD_Dis, (HC_Rum-HC_Dis)/HC_Dis
MRI_Phenotype = {};
for iCondition = 1: length(ConditionSet)
     load([TopOutputDir,'/ROISignals_', ConditionSet{iCondition}, '_lh.mat']);
     MRI_Phenotype{iCondition} = [network_FC{iCondition}, ROISignals];
end
MRI_Phenotype{5} = MRI_Phenotype{1} - MRI_Phenotype{2};
MRI_Phenotype{6} = MRI_Phenotype{3} - MRI_Phenotype{4};
MRI_Phenotype{7} = (MRI_Phenotype{1} - MRI_Phenotype{2})./MRI_Phenotype{2};
MRI_Phenotype{8} = (MRI_Phenotype{3} - MRI_Phenotype{4})./MRI_Phenotype{4};

MRI_label = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN core', ...
                            'DMN DMPFC', 'DMN MTL', 'SC', 'DC'};
condition_label = {'MDD_Rum', 'MDD_Dis', 'HC_Rum', 'HC_Dis', ...
                            'MDD_delta', 'HC_delta', ...
                            'MDD_relative_delta', 'HC_relative_delta'}; 
                        
%% do scale corr analyses
% MDD patients
stats_corr = [];
temp_MDD = [];
for i = [1,2,5,7]
    count = 0;
    for j = 1:size(MRI_Phenotype{1},2)
        for k = 1:size(scale_score,2)
            count = count + 1;
            y = MRI_Phenotype{i}(:,j);
            x = scale_score(Dx == 1, k);
            [r, p] = corr(x,y);
            eval(['stats_corr.',condition_label{i},'{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k}]']);
            eval(['stats_corr.',condition_label{i},'{count, 2} = r']);
             if p < 0.05
                eval(['stats_corr.',condition_label{i},'{count, 3} = p']);
                %store the significant condition ids
                temp_MDD(end+1,1) = i;
                temp_MDD(end,2) = j;
                temp_MDD(end,3) = k;
             end
        end
    end
end

% HC
temp_HC = [];
for i = [3,4,6,8]
    count = 0;
    for j = 1:size(MRI_Phenotype{1},2)
        for k = 1:size(scale_score,2)
            count = count + 1;
            y = MRI_Phenotype{i}(:,j);
            x = scale_score(Dx == 2, k);
            [r, p] = corr(x,y);
            eval(['stats_corr.',condition_label{i},'{count, 1} = [MRI_label{j}, ''-'', scale_name_set{k}]']);
            eval(['stats_corr.',condition_label{i},'{count, 2} = r']);
             if p < 0.05
                eval(['stats_corr.',condition_label{i},'{count, 3} = p']);
                %store the significant condition ids
                temp_HC(end+1,1) = i;
                temp_HC(end,2) = j;
                temp_HC(end,3) = k;
             end
        end
    end
end

%%  extract the raw numbers for plot
% numbers
% i condition, j MRI phenotype, k scale number
% for ii = 1:length(temp_MDD)
i = 2; j = 3; k = 2;
y = MRI_Phenotype{i}(:,j);
x = scale_score(Dx == 1, k);
output_matrix_MDD = [];
output_matrix_MDD = [x, y];
output_matrix_MDD = [output_matrix_MDD, ones(length(output_matrix_MDD),1)];

i = 4;
x = []; y = [];
y = MRI_Phenotype{i}(:,j);
x = scale_score(Dx == 2, k);
output_matrix_HC = [];
output_matrix_HC = [x, y];
output_matrix_HC = [output_matrix_HC, ones(length(output_matrix_HC),1).*2];

output_matrix = [];
output_matrix = [output_matrix_MDD; output_matrix_HC];

save([work_dir,'/DAN_Reflection_corr_V3.txt'], 'output_matrix', '-ASCII', '-tabs');
% end

%% do mediation analysis
addpath(genpath('/mnt/Data2/RfMRILab/Chenxiao/CX_software/MediationToolbox-master'));
addpath(genpath('/mnt/Data2/RfMRILab/Chenxiao/CX_software/CanlabCore-master'));

%% do explorative median analysis, warning: very time consuming!
%MDD
iCount = 0;
for iCondition = [1,2, 5, 7]
    for iX = 1:size(scale_score,2)
            for iM = 1:size(MRI_Phenotype{1},2)
                for iY = 1:size(scale_score,2)
                    if iX == iY
                        continue;
                    else
                         X = scale_score(Dx == 1, iX);
                         M = MRI_Phenotype{iCondition}(:,iM);
                         Y = scale_score(Dx == 1, iY);
                         [~, stats] = mediation(X, Y, M, 'verbose', 'boot', 'bootsamples', 10000);
                         if stats.p(1) <0.05 && stats.p(2) < 0.05 && stats.p(5) < 0.05
                            iCount = iCount + 1;
                            output_median.label{iCount, 1} = scale_name_set{iX};
                            output_median.label{iCount, 2} = MRI_label{iM};
                            output_median.label{iCount, 3} = scale_name_set{iY};
                            output_median.label{iCount, 4} = condition_label{iCondition};
                            output_median.stats_names = stats.names;
                            output_median.p(iCount,:) = stats.p;
                         end
                    end
                end
            end
    end
end

%% get the significant results
% MDD
X = []; M = []; Y = [];
X = scale_score(Dx == 1, 2);
M = MRI_Phenotype{2}(:,3);
Y = scale_score(Dx == 1, 8);

[paths, stats] = mediation(X, Y, M, 'plots', 'verbose', 'boot', 'bootsamples', 10000);

% HC
X = scale_score(Dx == 2, 2);
M = MRI_Phenotype{4}(:,3);
Y = scale_score(Dx == 2, 8);

[paths, stats] = mediation(X, Y, M, 'plots', 'verbose', 'boot', 'bootsamples', 10000);