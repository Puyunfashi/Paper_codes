% to conduct stats on regional metrics of Suzhou data
% Written by ChenXiao
% 20220331

%% initialization
clear; clc;
work_dir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/Regional_Metrics_HeadMotion_Updated';
if ~exist(work_dir, 'dir'); mkdir(work_dir); end
data_dir = '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Preprocessing';
TemplateDir = '/mnt/Data2/RfMRILab/Chenxiao/CX_software/DPABI-master/DPABISurf/SurfTemplates';

HemisphereNameSet = {'lh','rh'};
HemisphereSet = {'LH','RH'};
MeasureSet = {'ALFF','fALFF', 'ReHo', 'DegreeCentrality'};
% MeasureSet = {'DegreeCentrality'};
SuffixSet = {'','','','_Bilateral_PositiveWeightedSumBrain'};
% SuffixSet = {'_Bilateral_PositiveWeightedSumBrain'};
PipesuffixSet = {'_FunSurfWC','_FunSurfWC','_FunSurfWCF','_FunSurfWCF'};
% PipesuffixSet = {'_FunSurfWCF'};

%% select participants
demographic_info = readtable( ...
    '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Analysis/Demographic_Info_V3.txt','Delimiter','\t');

sub_list = demographic_info.Serial_Number;
Dx =  demographic_info.Dx;
Age = demographic_info.Age;
Sex = demographic_info.Sex;
Sex = Sex -1;
Edu = demographic_info.Edu;

%% get the scale scores
scale_score = [demographic_info.HAMD_Sum, demographic_info.HAMA_Sum, ...
                        demographic_info.BDI, demographic_info.Rumination, ...
                        demographic_info.Brooding, demographic_info.Reflection, ...
                        demographic_info.Childhood_Trauma];
scale_name_set = {'HAMD','HAMA','BDI','Rumination','Brooding','Reflection','Childhood_Trauma'};

%% some subjects do not have rum and dis data
% read in subjects who have rum/dis data
sub_list_rum_dis = importdata( ...
    '/mnt/Data2/RfMRILab/Chenxiao/Suzhou_Rumination/Preprocessing/Task/sub_list_task_merged20220118.txt');
% exclude those who don't have rum and dis data
WantedSubMatrix = ismember(sub_list, sub_list_rum_dis);
WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);

% scale_score = scale_score(WantedSubIndex,:);

%% rule out subjects with missing data, xx = 999
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

% scale_score = scale_score(WantedSubIndex,:);

%% Extract head motion
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

% scale_score = scale_score(WantedSubIndex,:);

%% rule out additional subjects
% Sub026 and Sub051's phase coding directions are different
% Sub005 is left handed
WantedSubMatrix = ones(length(sub_list),1);
WantedSubMatrix(strcmp(sub_list,  'Sub026')) = 0;
WantedSubMatrix(strcmp(sub_list,  'Sub051')) = 0;
WantedSubMatrix(strcmp(sub_list,  'Sub005')) = 0;

% % edu are large, exclude to balance
% WantedSubMatrix(strcmp(sub_list,  'Sub151')) = 0;
% WantedSubMatrix(strcmp(sub_list,  'Sub139')) = 0;
% WantedSubMatrix(strcmp(sub_list,  'Sub123')) = 0;
% WantedSubMatrix(strcmp(sub_list,  'Sub122')) = 0;
% 
WantedSubIndex = find(WantedSubMatrix);
sub_list = sub_list(WantedSubIndex);
Dx = Dx(WantedSubIndex);
Age = Age(WantedSubIndex);
Sex = Sex(WantedSubIndex);
Edu = Edu(WantedSubIndex);
HeadMotion = HeadMotion(WantedSubIndex,:);

% scale_score = scale_score(WantedSubIndex,:);

%% check whether demographical info is balanced
% test the balance of two groups' age, sex, edu
% age
[~,p,~,stats] = ttest2(Age(Dx == 1), Age(Dx == 2));
stats_dem(1,1) = stats.tstat;
stats_dem(1,2) = p;
% Edu
[~,p,~,stats] = ttest2(Edu(Dx == 1), Edu(Dx == 2));
stats_dem(2,1) = stats.tstat;
stats_dem(2,2) = p;
% sex
[~,chi2,p,~] = crosstab(Dx,Sex);
stats_dem(3,1) = chi2;
stats_dem(3,2) = p;
% head motion rumination
[~,p,~,stats] = ttest2(HeadMotion(Dx == 1,1), HeadMotion(Dx == 2,1));
stats_dem(4,1) = stats.tstat;
stats_dem(4,2) = p;
% head motion distraction
[~,p,~,stats] = ttest2(HeadMotion(Dx == 1,2), HeadMotion(Dx == 2,2));
stats_dem(5,1) = stats.tstat;
stats_dem(5,2) = p;

full_demo_data = [Age,Edu,HeadMotion];
length(find(Dx == 1))
length(find(Dx == 2))
sum(Sex(Dx == 1))
sum(Sex(Dx == 1))/length(find(Dx == 1))
sum(Sex(Dx == 2))
sum(Sex(Dx == 2))/length(find(Dx == 2))

demog_MDD(1,:) = mean(full_demo_data(Dx == 1,:), 1);
demog_MDD(2,:) = std(full_demo_data(Dx == 1,:));
demog_HC(1,:) = mean(full_demo_data(Dx == 2,:), 1);
demog_HC(2,:) = std(full_demo_data(Dx == 2,:));

%%  Arrange Files
%% Organize data, could be skipped
% Rumination
% MDD
for i = 1:length(find(Dx == 1))
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                              MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                              '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/rum/MDD_FunSurf',HemisphereSet{iHem},'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
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
                                 MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                                 '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/rum/HC_FunSurf',HemisphereSet{iHem},'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

% Distraction
%MDD
for i = 1:length(find(Dx == 1))
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/S2_ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                             MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                              '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/dis/MDD_FunSurf',HemisphereSet{iHem},'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%HC
for i = length(find(Dx == 1))+1:length(sub_list)
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/S2_ResultsS/FunSurf',HemisphereSet{iHem}, ...
                                    '/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                                    '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/dis/HC_FunSurf',HemisphereSet{iHem},'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%% do mixed effect analysis
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

% left hemisphere
PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_lh_white.surf.gii'];
PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_lh_white_avg.area.gii'];
MaskFile =  [TemplateDir,'/fsaverage5_lh_cortex.label.gii'];
for iMeasure = 1:length(MeasureSet)
        DependentDir{1,1} = [work_dir,'/rum/MDD_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{2,1} = [work_dir,'/dis/MDD_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{3,1} = [work_dir,'/rum/HC_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{4,1} = [work_dir,'/dis/HC_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        OutputDir = [work_dir,'/stats'];
        if ~exist(OutputDir); mkdir(OutputDir); end
        OutputName = [OutputDir,'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure},'_lh'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
end

% right hemisphere
PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_rh_white.surf.gii'];
PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_rh_white_avg.area.gii'];
MaskFile =  [TemplateDir,'/fsaverage5_rh_cortex.label.gii'];
for iMeasure = 1:length(MeasureSet)
        DependentDir{1,1} = [work_dir,'/rum/MDD_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{2,1} = [work_dir,'/dis/MDD_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{3,1} = [work_dir,'/rum/HC_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{4,1} = [work_dir,'/dis/HC_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        OutputDir = [work_dir,'/stats'];
        if ~exist(OutputDir); mkdir(OutputDir); end
        OutputName = [OutputDir,'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure},'_rh'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
end

%% extract significant cluster
% I found DC is significantly interacted in patients and HCs and the values are extracted for plotting 
Group_Set = {'MDD', 'HC'};
Condition_Set = {'rum', 'dis'};
counter = 0;
output_matrix = [];
for ii = 1:length(Group_Set)
    for jj = 1:length(Condition_Set)
        counter = counter + 1;
        ROISignals = [];
        ROIDef{1} = [work_dir, '/DC_interaction_Cluster_Mask.gii'];
        AllVolume = [work_dir, '/', Condition_Set{jj}, '/', Group_Set{ii}, '_FunSurfLH/DegreeCentrality_FunSurfWCF'];
        OutputName = [work_dir, '/DC_Interaction_', Group_Set{ii}, '_', Condition_Set{jj}, '_LH'];
        [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, '', 0);
        ROISignals = [ROISignals, ones(length(ROISignals),1).*counter];
        output_matrix = [output_matrix; ROISignals];
    end
end

save([work_dir,'/DC_Interaction.txt'], 'output_matrix', '-ASCII', '-tabs');

%% do post hoc comparisons
output_matrix = importdata([work_dir,'/DC_Interaction.txt']);
[~,p,~,stats] = ttest(output_matrix(output_matrix(:,2) == 1), output_matrix(output_matrix(:,2) == 2));
posthoc.mddrum_vs_mdddis(1) = stats.tstat;
posthoc.mddrum_vs_mdddis(2) = p;

[~,p,~,stats] = ttest(output_matrix(output_matrix(:,2) == 3), output_matrix(output_matrix(:,2) == 4));
posthoc.hcrum_vs_hcdis(1) = stats.tstat;
posthoc.hcrum_vs_hcdis(2) = p;

[~,p,~,stats] = ttest2(output_matrix(output_matrix(:,2) == 1), output_matrix(output_matrix(:,2) == 3));
posthoc.mddrum_vs_hcrum(1) = stats.tstat;
posthoc.mddrum_vs_hcrum(2) = p;

[~,p,~,stats] = ttest2(output_matrix(output_matrix(:,2) == 2), output_matrix(output_matrix(:,2) == 4));
posthoc.mdddis_vs_hcdis(1) = stats.tstat;
posthoc.mdddis_vs_hcdis(2) = p;


%% do paired t tests
TopOutput_Dir = [work_dir, '/stats/paired_t'];
if ~exist(TopOutput_Dir, 'dir'); mkdir(TopOutput_Dir); end

PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

% MDD
OtherCovariates = {};
OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);

for iMeasure = 1:length(MeasureSet)
    for iHem = 1:2
        PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white.surf.gii'];
        PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white_avg.area.gii'];
        MaskFile =  [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_cortex.label.gii'];
        DependentDirs{1,1} = [work_dir,'/rum/MDD_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDirs{2,1} = [work_dir,'/dis/MDD_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        OutputName = [TopOutput_Dir,'/PT_MDD_',HemisphereNameSet{iHem}];
        [TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,PALMSettings);
    end
end

% HC
OtherCovariates = {};
OtherCovariates{1,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{2,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

for iMeasure = 1:length(MeasureSet)
    for iHem = 1:2
        PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white.surf.gii'];
        PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_white_avg.area.gii'];
        MaskFile =  [TemplateDir,'/fsaverage5_',HemisphereNameSet{iHem},'_cortex.label.gii'];
        DependentDirs{1,1} = [work_dir,'/rum/HC_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDirs{2,1} = [work_dir,'/dis/HC_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        OutputName = [TopOutput_Dir,'/PT_HC_',HemisphereNameSet{iHem}];
        [TTestPaired_T,Header] = y_TTestPaired_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,PALMSettings);
    end
end