%ExtractImarisData
% Code to Retrieve and Colosolidate Morphology Data from Imaris CV files
% edited starting 6/11/24 Ashley (AMP) and Emma (ES)

%% Notes
% this code should iterate through a folder of subject files (one per subject) and populate a new data structure. 
% 
% Imaris data should be acquired by researchers blind to experimental conditions.
% The new data structure can be passed to the unblinding code within this script by turning on a line.
% 
% The script should then go through ANOVA tests,
% post hoc comparisons between groups if there are two or more groups, and graphs made with violin plot
% be sure to specify whether data is raw, or averaged or totaled by cell,
% as there is separate code to create averaged data

% YOU MUST SPECIFY THE REGION AND MEASURE OF INTEREST AS WELL AS IF DATA
% ARE BLINDED WITH DIFFERENT EXPERIMENTAL NAMES via using a Key

% If you have other variables, such as stress, or genotype, etc. you must
% add them as new columns and append the empty dummy as well as new tables 

%% Install and Operating Notes
% make sure to install Violin Plot https://github.com/bastibe/Violinplot-Matlab
% make sure Paths are updated
% create an experimental Key with names {'Subject' 'Condition' 'Region' 'Sample' };
% modify as need for things like genotype, sex, and blinded names if aplicable etc.
% 'Sample' or slice is necessary if there are more than one sample per animal per
% region (e.g. two slices were imaged of hippocampus per mouse and quantified in Imaris)

% clear everything that was open
clc;
clear;
close all;
disp ('if something goes wrong from here, its likely a path problem, a cannot find your data issue');
disp ('or because Matlab is stupid.');
disp(' ');

prompt_measures = {'Enter the location for the data source;'};
cellopts = {'Astrocytes', 'Microglia'};
[indx,tf]=listdlg('ListString', cellopts); % dialog box for options
celltype = string(cellopts(indx)); % get list number of option

% Set Up Paths - Ashley's Paths
% % directory locations as follows
addpath(pwd); % set principle working directory to the experiment folder
basepath = '\\ad.gatech.edu\bme\labs\singer';% %make a path to the cell morphology folder
logpathname = ([basepath, '\Ashley\Projects\Tina_Thy_Study\Morphology\']);
addpath(logpathname);
% 
% %make a directory of the experiment folder too
experiment_folder = ([basepath, '\Ashley\Projects\Tina_Thy_Study\Morphology']);
experiment_dir=dir(experiment_folder);
% 
% %location of code
codepathname = ([basepath, '\Ashley\Projects\Tina_Thy_Study\Morphology\Code\']);
addpath(codepathname);
% %location of data including all folders and subfolders for each animal
% datapathname = ([basepath, '\Tina\Original Thy1 images obf\mybox-selected\Statistical Data\Subjects\']);
% datapathname = ([basepath, '\Tina\For AMP_convex and filament_June2024_exported files\']);
datapathname = ([basepath, '\Tina\Original Thy1 images obf\Updated male microglia renderings_May 2025\Microglia\Filaments\PFC\']);

% %add violin plot folder to path
violinplot_location=([basepath, '\Ashley\Violinplot-Matlab-master']);
addpath(violinplot_location);

sigstar_location = ([basepath, '\Ashley\Projects\Tina_Thy_Study\Morphology\Code\github_repo']);
addpath(sigstar_location)

%% Save to paths
% % location of path to save new consolidated data to
if not(isfolder([logpathname 'Data\']))
    mkdir([logpathname 'Data\'])
end
newdatapath = ([logpathname 'Data\']);
addpath(newdatapath);

% set path to save figures to if they don't already exist
if not(isfolder([logpathname 'Graphs\']))
    mkdir([logpathname 'Graphs\'])
end
figpathname = ([logpathname 'Graphs\']);
addpath(figpathname);

%% edit below for region based on the data (ex HPC or VC or PFC)
%****************** specify REGION to analyze here **** run one measure at a time ***********
region = 'mPFC';

% datapathname=([datapathname region '\']);

disp('paths are set up');

% removing warning for variable names in tables 
id = 'MATLAB:table:ModifiedAndSavedVarnames';
warning('off',id);

%% Most of what is below will not need to be edited unless for your experimental needs
% or if you have other measurs of interest or difficult paths to edit

% you must edit the violin plots at the end of the script !!!

%% Load reference file and a new data structure to populate
% set up an experimental key and load it here
% key should have names {'Subject' 'Condition' 'Region' 'Sample' };

% expfiles = experiment_dir(~([experiment_dir.isdir])) % find all files not folders
% count=1;
% while (count<=length(expfiles))
%     expfilename = expfiles(count).name
%     foundKey = strfind(expfilename, 'Key') %look for file with "key" in name
%     if ~isempty(foundKey) % if the rsults are not empty, it has found the keyfile 
%         key_doc = expfilename
%     end
%     count = count + 1;
% end
key_doc = fullfile(pwd,'ChronicFlickerStress_MouseKey_TF_updated.xlsx'); %for this line of code, get out of Code and get back into Morphology folder where excel sheet is found
mousekey = readtable(fullfile(key_doc));
missingmice=ismissing(mousekey.ImageSubjectID);
mousekey(missingmice,:)= []; % remove mice from key that do not have image files
number_of_samples = height(mousekey); % how many mice in key file
disp('a key has been found');

% set up a new table to populate with the experimental data
% originfo = array2table(NaN(1,6)); % create a new info table
% originfo.Properties.VariableNames = {'Subject' 'Condition' 'Region' 'FilamentID' 'VariableofInterest' 'Sex'};

% use this line below if you have multiple samples for each mouse
originfo = array2table(NaN(1,10));
originfo.Properties.VariableNames = {'Subject' 'Condition' 'Region' 'Sample' 'FilamentID' 'VariableofInterest' 'Sex' 'State' 'Flicker' 'Overview'};

%% morphology specifics to look for in the Imaris file names from detailed statistics

% Soma Measures
soma_number = 'tbd'; % should be num_cells
soma_volume = 'Soma_Volume'; % soma or cell body volume um3 
dist_to_nearest_neighbor = 'Distance_To_Nearest_Neighbour';
nearest_neighbor_3 = 'Average_Distance_To_3_Nearest_Neighbours';
nearest_neighbor_5 = 'Average_Distance_To_5_Nearest_Neighbours';
nearest_neighbor_9 = 'Average_Distance_To_9_Nearest_Neighbours';

% surface & process measures
total_process_volume = 'Filament_Volume_(sum)'; % total process volume per microglia - the soma (no soma here)
total_process_length = 'Filament_Length_(sum)'; % total length per cell
branching_depth = 'Filament_Full_Branch_Depth';
branch_points = 'Filament_No._Segment_Branch_Pts'; % number of nodes
% Filament_No._Segment_Branch_Pts.csv
number_of_branches = 'Filament_No._Dendrite_Branch_Pts'; % number of branches resulting from nodes

% segment measures (node to node segments and not the entire process)
average_segment_length = 'Dendrite_Length';
average_segment_volume = 'Dendrite_Volume';
average_segment_level = 'Dendrite_Branch_Level';

% convex hull volume
% convexhull_volume = 'convex_Volume'; % convexHull_volume
convexhull_volume = 'convex_Detailed'; % convexHull_volume pulled from each row

% arborization is a calculated measure
% arborization = dendrite volume (sum) / filament No. dendrite branches (sum)

% colocalization unique to Tina's Thy study output
neuron_colocalization = 'interaction_Average';
microglia_neuron_colocalization = 'interaction_Average';
% row 46, column I = Overlapped Volume Ratio to Surfaces (neuron) Sum
% row 47, column I = Overlapped Volume Ratio to Surfaces (microglia neuron contact) Sum

disp('variables of interest are set');

%% ****************** specify string to find here **** run one measure at a time ***********
stringToBeFound = total_process_length;  %convexhull_volume

% below specific to Tina Thy study

coloc_calc= strcmp(stringToBeFound,'interaction_Average');
if coloc_calc == 1
    datapathname = ([basepath, '\Tina\For AMP_neuron interaction_June2024\']);
    % colocalization needs to be a string input for cell type
    prompt = "Colocalization for microglia-neurons [m] or just neurons [n]? ";
    txt = input(prompt,"s");
    coloc_input = strcmp(txt,"m");

    if coloc_input == 1
        coloc_target = 1; % m = 1 for microglia
    else
        coloc_target = 2; % otherwise n = 2 for neurons
        txt = 'Y';
    end
else
    coloc_target = 0;
end

% filament_calc = contains(stringToBeFound, 'Filament');
% if filament_calc == 1
%     datapathname = ([basepath, '\Tina\For AMP_filaments_June2024\']);
% end
% 
% soma_calc = contains(stringToBeFound, 'Soma');
% if soma_calc == 1
%     datapathname = ([basepath, '\Tina\For AMP_filaments_June2024\']);
% end

%%make an array for stringtobefound
%%arrayToBeFound = [];
disp(['now calculating ' stringToBeFound ' for each sample in the ' region ' using means and maxes for each Filament ID...']); % sanity check
% % uncomment above for sanity check otherwise leave off

%% read in folders and begin iteration through folder
filesAndfolders = dir(datapathname); % returns all files and folders in directory
numoffolders = length(filesAndfolders); % the number of folders in the directory
b = 1:size(filesAndfolders);
i=3; % set a folder to start, first 3 folders are filler '..' and should be skipped for this process

averages_table=table; % create table to update each iteration of averages
updatedinfo = originfo; %create table to update each iteration of raw values

disp('images processed after the 4th image ')
while(i<=numoffolders)
    curdir = filesAndfolders(i).name;% Store the name of the folder
    curdir = lower(curdir);% Make it case insensitive
%     subjectID=char(curdir);
    subjectID=extractBefore(curdir, "_");
%     subjectID=extractBetween(subjectID, "_", "_slice", 'boundaries', 'exclusive');
    subjectID=string(subjectID);
    beginningoftile = strfind(curdir, 'tile');
    sampleID = curdir(beginningoftile + 5); % get the two values after
    sampleID=string(sampleID);
%     sampleID = extractBetween(curdir, "tile_", "_", 'Boundaries','exclusive'); %%%come back to and search for tile_ and sample id will be # after _
    % add sex in here too
    % sampleID = string(sampleID);
    %disp(['for slice ' sampleID]); % sanity check comment
    overview = extractBetween(curdir,'overview ', '-', 'boundaries', 'exclusive');
    overview = string(overview{1,1});
    mouseinfo = ([datapathname curdir]); %create a path to the folder
    mousedir = dir(mouseinfo); % make it a directory
    num_dir=numel(mousedir); % count number of files in the directory
    all_filesandfoldersindir = dir(fullfile(mouseinfo,'/**/*.*'));
    filesindir = all_filesandfoldersindir(~([all_filesandfoldersindir.isdir])); % not folders
    numberoffiles = length(filesindir); % the number of files in the directory
    
    a=1; %create an interation
    while(a<=numberoffiles) %start the while loop through files within the directory
        filename = filesindir(a).name;  % Store the name of the file
        found = strfind(filename, stringToBeFound); % look for specific variable in file name listed above
        if ~isempty(found) % if the results are not empty, it has found the file name
            foundString = strcat('Found in file ------', filename);
%           disp(foundString); % sanity check
            thismouse = [datapathname curdir];      
            % filename = filename(~ismember(filename, ']['));%remove brackets from filename
            % disp(filename)
            
            data = readtable(fullfile(thismouse,filename)); % read in excel with subject data
            
            % *********** specific to Tina Thy Study  *******
            
            coloc_calc = contains(stringToBeFound, 'interaction_Average');
            
            if coloc_calc == 1
                cvdata = readtable(fullfile(thismouse,filename),'ReadVariableNames', true ); % read in excel with subject data
                tempcoloc_idx = contains(cvdata.Variable, 'Overlapped Volume Ratio to Surfaces');
                colocinfo= cvdata(tempcoloc_idx, :);
 
                if coloc_target == 1 % microglia-neuron surface ratio contact sum
                    glia_neuron_coloc_calc = contains(colocinfo.Surfaces, 'microglia'); % sanity check
                    if sum(glia_neuron_coloc_calc) == 1
                        %                         data = readtable(fullfile(thismouse,filename),'ReadVariableNames', true ); % read in excel with subject data
                        %                         tempcoloc_idx = contains(data.Variable, 'Overlapped Volume Ratio to Surfaces');
                        %                         colocinfo= data(tempcoloc_idx, :);
                        midx = contains(colocinfo.Surfaces, 'microglia');
                        data = colocinfo(midx,:); % microglia-neuron data is second row
                        data.ID = "microglia-neuron contact";
                        data.(1) = data.Sum;
                        neuron_coloc_calc = 0;
                    else
                        data = table(1,1);
                        data.ID = "microglia-neuron contact";
                        data.(1) = 0;
                    end
                elseif coloc_target == 2
%                     neuron_coloc_calc = contains(colocinfo.Surfaces, 'neuron'); % sanity check
                    neuron_coloc_calc = strncmp(colocinfo.Surfaces,'neuron',6);
                % row 46, column I = Overlapped Volume Ratio to Surfaces (neuron) Sum
                    if sum(neuron_coloc_calc) == 1
%                     tempcoloc_idx = contains(data.Variable, 'Overlapped Volume Ratio to Surfaces');
%                     colocinfo= data(tempcoloc_idx, :);
                        nidx=strncmp(colocinfo.Surfaces,'neuron',6);
                        data = colocinfo(nidx,:);
                        data.ID = "neuron contact";
                        data.(1) = data.Sum;
                        glia_neuron_coloc_calc = 0;
                    else
                        data=table(1,1);
                        data.ID = "neuron contact";
                        data.(1) = 0;
                    end
                else
                    glia_neuron_coloc_calc = 0;
                    neuron_coloc_calc = 0;
                end
            end
            % turn on below to replace outliers within sample data
            %data.(1) = filloutliers(data.(1),NaN); %replace outliers from median with NaNs
            % num_cells = height(data); % number of filament rows = cell count
            num_cells = (height(data)) ; % number of filament rows = cell count 
            dummy_mouse = repmat(subjectID,num_cells,1); % populate mouse name
            dummy_condition = NaN(num_cells,1); % populate cell data
            dummy_region = repmat(string(region),num_cells,1); % populate cell data
            dummy_sample = repmat(sampleID, num_cells,1); % populate sample number data
            dummy_overview = repmat(overview, num_cells,1); % populate sample number data
            dummy_sex = NaN(num_cells,1);%
            dummy_state = NaN(num_cells,1); % populate genotype data
            dummy_flicker = NaN(num_cells,1); % populate flicker frequency data
            newinfo = table(dummy_mouse, dummy_condition, dummy_region, dummy_sample, data.ID, data.(1), dummy_sex, dummy_state, dummy_flicker, dummy_overview);
            newinfo.Properties.VariableNames = {'Subject' 'Condition' 'Region' 'Sample' 'FilamentID' 'VariableofInterest' 'Sex' 'State' 'Flicker' 'Overview'}; % filament ID = cell ID
%             disp 'obtaining sample info'; % sanity check

%% get average data for soma volume
            avgcalc= strcmp(stringToBeFound,'Soma_Volume');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Filament_maxvols = splitapply(@max,data.SomaVolume,G); % find maximum dendrite length for each filament ID
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'SomaVolume'});
                meanbycell.MaxSomaVolume = Filament_maxvols;
                h = height(meanbycell);
                temp = {subjectID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessVolume' 'MaxProcessVolume'};
                averages_table=[averages_table; avginfo];

                updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
%                 save(updated_avgfile,'averages_table');
%                 disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
            end
         %% get average data for convex hull area
            avgcalc= strcmp(stringToBeFound,'convexHull_Area');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Filament_maxvols = splitapply(@max,data.convexHullArea,G); % find maximum dendrite length for each filament ID
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'convexHullArea'});
                meanbycell.convexHullArea = Filament_maxvols;
                h = height(meanbycell);
                temp = {subjectID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessVolume' 'MaxProcessVolume'};
                averages_table=[averages_table; avginfo];

                updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
%                 disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
            end
          %% get average data for convex hull volume
    %         avgcalc= strcmp(stringToBeFound,'convex_Detailed');
    %         if avgcalc == 1
    %             [G,maxlen] = findgroups(data.FilamentID);
    %             ConvexHull_maxvols = splitapply(@max,data.ConvexHullVolume,G); % find maximum dendrite length for each filament ID
    %             omean = @(x) mean (x,'omitnan');
    %             meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'convex_Detailed'});
    %             meanbycell.MaxConvexHullVolume = ConvexHull_maxvols;
    %             h = height(meanbycell);
    %             temp = {subjectID, region}; % create a table of new info
    %             avginfo = repmat(temp, h,1);
    %             avginfo = cell2table(avginfo);
    %             avginfo = [avginfo meanbycell];
    % % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
    %             avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessVolume' 'MaxProcessVolume'};
    %             averages_table=[averages_table; avginfo];
    % 
    %             updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
    %             save(updated_avgfile,'averages_table');
    %             disp(['saved averages data ' updated_avgfile]);
    %             save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
    %         end
            %% use below to get means by groups for dendrite length
            avgcalc= strcmp(stringToBeFound,'Dendrite_Length');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Dendrite_maxlens = splitapply(@max,data.DendriteLength,G); % find maximum dendrite length for each filament ID                
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables', 'Sample''FilamentID','InputVariables',{'DendriteLength'});
                meanbycell.MaxDendriteLength = Dendrite_maxlens;
                h = height(meanbycell);
                temp = {subjectID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};
                averages_table=[averages_table; avginfo];
                updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
                disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
            end
     %% use below to get means by groups for filament area
            avgcalc= strcmp(stringToBeFound,'Filament_Area_(sum)');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Filament_maxarea = splitapply(@max,data.FilamentArea,G); % find maximum dendrite length for each filament ID                
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables', 'Sample''FilamentID','InputVariables',{'DendriteLength'});
                meanbycell.MaxFilamentArea = Filament_maxarea;
                h = height(meanbycell);
                temp = {subjectID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};
                averages_table=[averages_table; avginfo];
                updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
                disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
            end
            %% get average data for segment spine density
            avgcalc= strcmp(stringToBeFound,'Segment_Spine_Density');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Segment_maxdens = splitapply(@max,data.SegmentSpineDensity,G); % find maximum dendrite length for each filament ID
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'SegmentSpineDensity'});
                meanbycell.MaxSegmentSpineDensity = Segment_maxdens;
                h = height(meanbycell);
                temp = {subjectID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessVolume' 'MaxProcessVolume'};
                averages_table=[averages_table; avginfo];

                updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
                disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
            end
% %             %% use below to get means by groups for Soma Area
% %             avgcalc= strcmp(stringToBeFound,'Soma_Area');
% %             if avgcalc == 1
% %                 [G,maxlen] = findgroups(data.FilamentID);
% %                 Soma_maxarea = splitapply(@max,data.SomaArea,G); % find maximum soma area for each filament ID                
% %                 omean = @(x) mean (x,'omitnan');
% %                 meanbycell = varfun(omean,data,'GroupingVariables', 'Sample''FilamentID','InputVariables',{'SomaArea'});
% %                 meanbycell.MaxSomaArea = Soma_maxarea;
% %                 h = height(meanbycell);
% %                 temp = {subjectID, region}; % create a table of new info
% %                 avginfo = repmat(temp, h,1);
% %                 avginfo = cell2table(avginfo);
% %                 avginfo = [avginfo meanbycell];
% %     % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
% %                 avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessLength' 'MaxProcessLength'};
% %                 averages_table=[averages_table; avginfo];
% %                 updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
% %                 save(updated_avgfile,'averages_table');
% %                 disp(['saved averages data ' updated_avgfile]);
% %                 save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
% %             end
            %% get average data for dendrite volume
            avgcalc= strcmp(stringToBeFound,'Dendrite_Volume');
            if avgcalc == 1
                [G,maxlen] = findgroups(data.FilamentID);
                Filament_maxvols = splitapply(@max,data.DendriteVolume,G); % find maximum dendrite length for each filament ID
                
                omean = @(x) mean (x,'omitnan');
                meanbycell = varfun(omean,data,'GroupingVariables','FilamentID','InputVariables',{'DendriteVolume'});
                meanbycell.MaxDendriteVolume = Filament_maxvols;
                h = height(meanbycell);
                temp = {subjectID, region}; % create a table of new info
                avginfo = repmat(temp, h,1);
                avginfo = cell2table(avginfo);
                avginfo = [avginfo meanbycell];
    % %             newinfo.Properties.VariableNames = {'Mouse' 'Region' 'Sample' 'FilamentID' 'NumberofSegments' 'MeanProcessLength'};
                avginfo.Properties.VariableNames = {'Subject' 'Region' 'FilamentID' 'NumberofSegments' 'MeanProcessVolume' 'MaxProcessVolume'};
                averages_table=[averages_table; avginfo];

                updated_avgfile=([logpathname 'Data/' region '_averages_' stringToBeFound '_appended_' date '.mat']);
                save(updated_avgfile,'averages_table');
                disp(['saved averages data ' updated_avgfile]);
                save([logpathname 'Data/' stringToBeFound '_averages_table.mat'],'averages_table');
            end
            updatedinfo = [updatedinfo; newinfo]; % append rows to table for each subset of mouse data

            %% use below to get means by cell for each subject if desired
%             avgdata=grpstats(newinfo,["Subject","FilamentID"],["mean","sem"],"DataVars",'VariableofInterest');
%             averages_table = [averages_table; avgdata];
%             disp 'getting cell based averages'; % sanity check
            break;
        end
        a=a+1; % iterate through additional files in folder if not found
    end
    i = i+1; % move onto the next mouse ID
    fprintf('%d ',i) % display iteration number sanity check
end        
    
%%        
variablename = strrep(stringToBeFound,'_', ' '); %repalce _ with ' ' in filename

% *** specific for Tina Thy study ***
if coloc_target == 1
    variablename = 'microglia neuron contact';
    neuron_coloc_calc = 0;
end
if coloc_target == 2
    variablename = 'neuron contact';
    glia_neuron_coloc_calc = 0;
end
% ***********

if height(updatedinfo) < 1
    error_msg = ' Something went wrong! No files were found.';
    error_warning = msgbox(error_msg);
    error(error_msg)
end

disp(' ');
disp(['creating a structure called updatedinfo for ' variablename]); % display sanity check for subject
disp(['for the region ' region]); %display sanity check for region
disp(' ');

updatedinfo(1,:) = []; % delete first blank row
updatedinfo=renamevars(updatedinfo,"VariableofInterest",variablename);

%% Pass data to unblinding code - if using researchers blind to conditions also copy unblinded names
% the output is "updated_unblindedinfo" below - updated 6/13/24 AMP to
% account for sex as a variable

disp('passing raw data to unblinding code using the key');  

updated_unblindedinfo=updatedinfo; % rename a new table with updated info for saving

keytable_height = height(mousekey); % how many mice in key file
updatedtable_height = height(updatedinfo); % how much data to parse
condition_placeholder=string(nan(updatedtable_height,1)); % dummy array of Nans to populate for condition
blind_name_placeholder = string(nan(updatedtable_height,1)); % dummy array of Nans to populate for real names
sampleID_placeholder = string(nan(updatedtable_height,1)); % dummy array of Nans for populate for sample IDs
sex_placeholder = string(nan(updatedtable_height,1)); % dummy array of Nans to populate for sex
flicker_placeholder = string(nan(updatedtable_height,1)); % dummy array of Nans to populate for flicker
state_placeholder = string(nan(updatedtable_height,1)); % dummy array of nans to populate for state

x = 1; % start iteration at 1 through the mouse key
while(x<=keytable_height)
    % keymouse = string(mousekey.Subject(x)); % for this mouse in the key table
    % TURN ON line below if the mouse has a different blinded name on Key
    % keymouse = string(mousekey.ImageSubjectID(x)); % for this mouse's Blinded Imaging ID in the key table
    keymouse = string(mousekey.BlindSubjectID(x)); %5 May 5 TINA MALE MICROGLIA
    % keymouse = string(mousekey.PullSubjectID(x)); % for this mouse's Blinded Experimental ID in the key table
    % keymouse = string(mousekey.BlindName(x)); % for this mouse's Blinded Experimental ID in the key table
    updatedinfo.Subject = string(updatedinfo.Subject);
    idx = ismember(lower(updatedinfo.Subject),lower(keymouse)); % find index of mousekey name that matches data sample
    num_matches = sum(idx); % number of matching rows from data
    matching_data_row_num=find(idx); % row index numbers from data that have mouse
    condition = string(mousekey.Condition(x)); % what the experimental group was
    condition_placeholder(matching_data_row_num) = condition; % fill in condition data per mouse
    blind_name_placeholder(matching_data_row_num) = keymouse; % fill in real names of subjects
    % sampleID_placeholder(matching_data_row_num) = string(mousekey.Sample(x));
    sampleID_placeholder(matching_data_row_num) = string(sampleID);
    sex= string(mousekey.Sex(x)); % what the sex was 
    sex_placeholder(matching_data_row_num) = sex; % fill in sex data per mouse
    flicker = string(mousekey.Flicker(x)); % what the flicker frequency was 
    state= string(mousekey.State(x)); % what the state was
    flicker_placeholder(matching_data_row_num) = flicker; % fill in real flicker value per subject
    state_placeholder(matching_data_row_num) = state; % fill in real state value per subject

    x = x+1; % increase the iteration
end

% update information on spreadsheet from placeholder info
updated_unblindedinfo.Condition=condition_placeholder; 
updated_unblindedinfo.ExperimentalID=blind_name_placeholder;
updated_unblindedinfo.Sex=sex_placeholder;
updated_unblindedinfo.Flicker=flicker_placeholder;
updated_unblindedinfo.State=state_placeholder;

updated_unblindedinfo(updated_unblindedinfo.(6) <= 0, : ) = []; % remove zero values for branch depth and number bifrucations

% updated_unblindedinfo=updatedinfo; % rename a new table with updated info for saving
disp('unblinding successful and data saved');

updatedfile=(strcat(logpathname, 'Data/', region, '_', variablename, '__appended_', date, '.mat'));
save(updatedfile,'updated_unblindedinfo'); % save the updated data structure under Data folder 
% maskchecked = updated_unblindedinfo.(6)==0;
% updated_unblindedinfo{maskchecked, 6}=NaN;
csvfile_update=(strcat(logpathname, 'Data/', region, '_', variablename, '__appended_', date, '.xls'));
writetable(updated_unblindedinfo,csvfile_update); % save info to excl file
% save([logpathname 'Data\updatedinfo.mat'],'updated_unblindedinfo');
% disp(['saved data structure ' updatedfile]);

%% Find Outliers - AMP edited 12/1/23
% disp('now checking for outliers within group and removing them');
% checked_data=table; %create empty table of checked data for outliers
% 
% conditionlist = unique(updated_unblindedinfo.Condition); % find unique groups (exp conditions)
% conditionlist=reshape(conditionlist,1,length(conditionlist)); % reshape for iteration
% conditionlist = rmmissing(conditionlist);
% number_of_conditions = length(conditionlist); % get number of exp groups
% % number_of_conditions = 6; %hard coded
% k = 1;
% while k<=length(conditionlist) % iterate through condition list
%     group=conditionlist(k); % specify group
%      conditiondata=ismember(updated_unblindedinfo.Condition,group); % find data for that group only
%      selected_cond_data=updated_unblindedinfo(conditiondata,:);  % create a table of that group's data
%      %%turn string into numerical value (cell2mat) and make sure its taking in the
%      %%right values
%      selected_cond_data.(8) = isoutlier(selected_cond_data.(6)); % create a column for the censor
%      selected_cond_data.Properties.VariableNames(8) = "OutlierCensor"; % label the censor
%      selected_cond_data.(6) = filloutliers(selected_cond_data.(6),NaN); %replace outliers from median with NaNs
%      checked_data=[checked_data; selected_cond_data]; %update table with outlier checked values
%      k=k+1;
% end
% 
% updated_unblindedinfo = checked_data; %repalce table with outliers removed 

% % special censor criteria for gial soma
% % for convex hull volume
% convex_vol_calc= strcmp(stringToBeFound,'convex_Detailed');
% if convex_vol_calc == 1
%     updated_unblindedinfo(updated_unblindedinfo.(6) <= 4.99, : ) = []; % remove soma volume less than 5um
% end
% 
% updatedfile=([logpathname 'Data/' region '_' variablename '__appended_' date '.mat']);
% save(updatedfile,'updated_unblindedinfo'); % save the updated data structure under Data folder 
% 
% csvfile_save = ([updatedfile '.xls']);
% writetable(updated_unblindedinfo,csvfile_save); % save info to excl file
% 
% % get averages per experimental condition
% conditionaverages=grpstats(updated_unblindedinfo,"Condition",["mean","sem"],"DataVars",variablename);
% conditionaveragesfile=([logpathname 'Data/' region '_' variablename '__Averages_appended_' date '.mat']);
% save(conditionaveragesfile,'conditionaverages'); % save the updated data structure under Data folder 

%% Clustering - EMS edited 6/16/2024 to account for Males and Females
% sex is assessed and graped separately for each measure of interest
% through an interation below to create separate tables

disp('');
disp('Separating out data based on sex');

maletable = array2table(NaN(1,8));
maletable.Properties.VariableNames = {'Subject' 'Condition' 'Region' 'Sample' 'FilamentID' 'VariableofInterest', 'Sex', 'Outlier Censor'};
femaletable = array2table(NaN(1,8));
femaletable.Properties.VariableNames = {'Subject' 'Condition' 'Region' 'Sample' 'FilamentID' 'VariableofInterest', 'Sex', 'Outlier Censor'};
malelist = [];
femalelist = [];
%create Male and Female tables

cellupdated = table2cell(updated_unblindedinfo);
for ij =1: height(updated_unblindedinfo)
    if cellupdated{ij,7} == 'M'
        malelist = [malelist, ij]; 
    else
        femalelist = [femalelist,ij];
    end
end
maletable = updated_unblindedinfo(malelist, :);
femaletable = updated_unblindedinfo(femalelist, :);

%% Statistics
% Check if the Distribution is normal nor not to determine parametric test
% use AD test to see if null is false (ouput = 1) or true (output = 0)
% if null is true (0), distribution is not normal and run rank sum test (if 2 groups or less)
% or KW test if more than 2 experimental groups
disp('');
disp('Beginning Statistics Now');
disp('testing data distributions for normality');


    %% Find Outliers ******************************************************** - AMP edited 8/14/24 to account for males and females
% disp('now checking for outliers within group and sex and removing them');
% checked_data=table; %create empty table of checked data for outliers
% updated_unblindedinfo_outliersremoved = table;
for i  = 1 % change to 1:2 if running for males and females
    if i == 1
        testgroup  = maletable;
        testgroupname = 'males';
    else
        testgroup = femaletable;
        testgroupname = 'females';
    end
% 
%     disp('now checking for outliers within group and removing them');
%     checked_data=table; %create empty table of checked data for outliers
% 
    conditionlist = unique(testgroup.Condition); % find unique groups (exp conditions)
    conditionlist = reshape(conditionlist,1,length(conditionlist)); % reshape for iteration
    conditionlist = rmmissing(conditionlist);
    number_of_conditions = length(conditionlist); % get number of exp groups
%     k = 1;
%     while k<=length(conditionlist) % iterate through condition list
%         group=conditionlist(k); % specify group
%         conditiondata=ismember(testgroup.Condition,group); % find data for that group only
%         selected_cond_data=testgroup(conditiondata,:);  % create a table of that group's data
%         %%turn string into numerical value (cell2mat) and make sure its taking in the
%         %%right values
%         selected_cond_data.(12) = isoutlier(selected_cond_data.(6), "median"); % create a column for the censor and specify method
%         selected_cond_data.Properties.VariableNames(12) = "OutlierCensor"; % label the censor
%         selected_cond_data.(6) = filloutliers(selected_cond_data.(6),NaN); %replace outliers from median with NaNs
%         checked_data=[checked_data; selected_cond_data]; %update table with outlier checked values
%         k=k+1;
%     end
%     updated_unblindedinfo_outliersremoved = [updated_unblindedinfo_outliersremoved; checked_data]; % add in outlier data column 12
% 
%     if i == 1
%         maletable = checked_data;
%         testgroup  = maletable;
%     else
%         femaletable = checked_data;
%         testgroup = femaletable;
%     end
%     % special censor criteria for gial soma
%     % for convex hull volume
%     convex_vol_calc= strcmp(stringToBeFound,'convex_Detailed');
%     if convex_vol_calc == 1
%         updated_unblindedinfo(updated_unblindedinfo.(6) <= 4.99, : ) = []; % remove soma volume less than 5um
%     end
%     updated_unblindedinfo = updated_unblindedinfo_outliersremoved; % replace table with outliers information
%     updatedfile=([logpathname 'Data/' region '_' variablename '__appended_outliersremoved_' date '.mat']);
%     save(updatedfile,'updated_unblindedinfo'); % save the updated data structure under Data folder 
% 
%     csvfile_save = ([updatedfile '.xls']);
%     writetable(updated_unblindedinfo,csvfile_save); % save info to excl file
% 
%     % get averages per experimental condition
%     conditionaverages=grpstats(updated_unblindedinfo,"Condition",["mean","sem"],"DataVars",variablename);
%     conditionaveragesfile=([logpathname 'Data/' region '_' variablename '__Averages_appended_' date '.mat']);
%     save(conditionaveragesfile,'conditionaverages'); % save the updated data structure under Data folder 


%% Iterate through male and female stats

    normal_group_test_results = array2table(NaN(1,2)); % create a table for test results
    normal_group_test_results.Properties.VariableNames = {'Group' 'normalitytestvalue'};
    normal_group_iter = 1; % create an iteration
    while normal_group_iter<=number_of_conditions
        test_condition = conditionlist(normal_group_iter); % find the condition from the list
        test_condition_data = testgroup.Condition == test_condition; % select that data
        normality_test_result = adtest(testgroup.(variablename));
        test_info = table(test_condition, normality_test_result);
        test_info.Properties.VariableNames = {'Group' 'normalitytestvalue'};
        normal_group_test_results=[normal_group_test_results; test_info];
        normal_group_iter = normal_group_iter + 1;
    end

    normal_group_test_results(1,:) = []; % delete first blank row
    sum_of_normality_tests = sum(normal_group_test_results.normalitytestvalue); % get the sum of the test scores
    % each group's score should = 0 if normal, = 1 if not normal
    % so the sum of scores should be greater than 0 if any group is not normal
    disp('Anderson-Darling Test for Normality performed for each group');
    disp('The AD result "h" is 1 or more if the test rejects the null hypothesis that data is from a population with a normal dist at the 5% significance level');
    norm_test_value = sprintf('The sum of the normality test values was %s', string(sum_of_normality_tests));
    disp(norm_test_value);
    test_cond_num = sprintf('The number of groups is %s ', string(number_of_conditions));
    disp(test_cond_num);
    % normality_val = logical(sum_of_normality_tests == number_of_conditions); % if they are equal (1) it is normal, if they are unequal (0) they are not normal
    normality_val = logical(sum_of_normality_tests > 0 ); % if they are zero it is normal, if they are greater than 0 they are not normal
    if normality_val == 0
        disp('The AD test null hypothesis was true, the distribution for groups was normal');
    else
        disp('The AD test null hypothesis was false, the distribution for one or more groups was not normal') ;
    end
    disp(' ');

    %clf;
    figure;
    % Statistical tests below depend on normality and number of conditions
    disp('determining appropriate statistics now...');
    disp('');
    disp(strcat( 'For ', testgroupname, ' data:'));
    if normality_val == 1 % groups were not normally  distributed based on test values, 1 = normally distributed
        % if sum of normality test scores is not equal to the number of groups, is not normal
        if number_of_conditions==2 % check if only 2 groups
            disp('Because the distribution for one or more groups was not normal, Wilcoxon Rank Sum Test was performed.'); % - not normal distribution test
            test1_condition = conditionlist(1); % find the 1st condition from the list
            test2_condition = conditionlist(2); % find the 2nd condition from the list
            test1_condition_data = testgroup.Condition == test1_condition; % select that data
            test2_condition_data = testgroup.Condition == test2_condition; % select that data
            [p,h,stats] = ranksum((testgroup.(variablename)(test1_condition_data)),(testgroup.(variablename)(test2_condition_data))); % run Wilcoxon Rank Sum Test
            wrfile = ([logpathname '\Stats\WRSum_Test_' variablename '_' testgroupname '_' date ]); % wr file save location
            save(wrfile,'stats'); % save wr test output
            text_multiflie = strcat(wrfile, '.txt');
            writestruct(stats, text_multiflie, "FileType", 'xml');
            % if p<0.05
            %     t_results = sprintf( 'For %s in the %s region, a Wilcoxon Rank-Sum test was run. %s was significantly different across groups (H (%s,%s)= %s, p = %s) for %s', celltype, region, variablename, string(t{2,3}), string(t{3,3}), string(t{2,6}), string(t{2,7}), testgroupname);
            %     disp(t_results);
            % else
            %     t_results = sprintf('For %s in the %s region, The results for %s were not significant for %s %s', celltype, region, variablename, genotype, testgroupname) ;
            %     disp(t_results);
            %     % results_text = (['Across conditions, there was a significant difference between groups (F (' string(t{2,3}) ',' string(t{3,3}) ') = ' string(t{2,6}) ' , p = ' string(t{2,7}) ';' ]);
            %     % The results are significant! The degree of freedom was %f the value of p is %f and the value of f was %f', t{2,3}, p,t{2,6}) ;
            % end
        else
            disp('Because the distribution for one or more groups was not normal, a KW test was performed.'); % - not normal distribution test for more than 2 groups
            % Run KW test Statistics for non normally distributed data with
            % more than 2 groups
            [p,r,stats] = kruskalwallis(testgroup.(6),testgroup.Condition); % KW test
            kwfile = strcat(logpathname, 'Stats\KW_Test_', variablename, '_', testgroupname, '_', date ); % kw file save location
            kwfile = string(kwfile);
            save(kwfile, 'r'); % save KW test output
            r = cell2table(r);
            text_multiflie = strcat(kwfile, '.txt');
            writetable(r, text_multiflie,'Delimiter','\t','WriteRowNames',true);

            if p<0.05
                t_results = sprintf( 'For %s in the %s region, a KW test was run. %s was significantly different across the rank totals for groups (H (%s)= %s, p = %s) for %s', celltype, region, variablename, string(r{2,3}), string(r{2,5}), string(r{2,6}), testgroupname);
                disp(t_results);
            else
                t_results = sprintf('For %s in the %s region, the results in for %s were not significant for %s %s', celltype, region, variablename, genotype, testgroupname) ;
                disp(t_results);
                % results_text = (['Across conditions, there was a significant difference between groups (F (' string(t{2,3}) ',' string(t{3,3}) ') = ' string(t{2,6}) ' , p = ' string(t{2,7}) ';' ]);
                % The results are significant! The degree of freedom was %f the value of p is %f and the value of f was %f', t{2,3}, p,t{2,6}) ;
            end

            % Run dunn's posthoc comparisons between groups
            [results,~,~,gnames] = multcompare(stats, 'CriticalValueType', 'dunn-sidak'); %returns posthoc paired comparisons test using correction for multiple comparisons
            dunn_posthocs= array2table(results,"VariableNames", ["Group","Control Group", "Lower Limit", "Difference", "upper Limit", "P-value"]);
            dunn_posthocs.("Group") = gnames(dunn_posthocs.("Group")); % identify the group names
            dunn_posthocs.("Control Group") = gnames(dunn_posthocs.("Control Group")); % specify the comparsion groups column 2
            multifile = ([logpathname '/Stats/KW_posthoc_comparisons_' variablename '_' region '_' testgroupname '_' date]);
            save (multifile, 'dunn_posthocs');
            text_multiflie = ([multifile '.txt']);
            writetable(dunn_posthocs, text_multiflie,'Delimiter','\t','WriteRowNames',true);
            posthocs_comps = dunn_posthocs;
        end
    
    elseif normality_val == 0 % normally distributed based on test values, 0 = normally distributed
        % else means that the sum of normality scores = the number of groups
        if number_of_conditions==2 % check if ONLY 2 groups
            disp('The distribution for all groups was normal, running a independent samples t test');
            test1_condition = conditionlist(1); % find the 1st condition from the list
            test2_condition = conditionlist(2); % find the 2nd condition from the list
            test1_condition_data = testgroup.Condition == test1_condition; % select that data
            test2_condition_data = testgroup.Condition == test2_condition; % select that data
            [h,p,~,stats] = ttest2((testgroup.(variablename)(test1_condition_data)),(testgroup.(variablename)(test2_condition_data)));
            tfile = ([logpathname '/Stats/T_Test_' variablename '_' testgroupname '_' date ]); % kw file save location
            save(tfile,'h'); % save t test output
            text_multiflie = ([tfile '.txt']);
            writetable(h, text_multiflie,'Delimiter','\t','WriteRowNames',true);
            if p<0.05
                t_results = sprintf( 'For %s in the %s region, a t-test was run. %s was significantly different across 40 Hz, 20 Hz and No Flicker groups (t (%s,%s)= %s, p = %s) for %s %s', celltype, region, variablename, string(t{2,3}), string(t{3,3}), string(t{2,6}), string(t{2,7}), genotype, testgroupname);
                disp(t_results);
            else
                t_results = sprintf('For %s in the %s region, the results for %s were not significant for %s %s', celltype, region, variablename, genotype, testgroupname) ;
                disp(t_results);
            end
        else
            disp('the distribution for all groups was normal, running an ANOVA');
            % Run Statistics - ANOVA if MORE than 2 groups
            [p,t,stats] = anovan(testgroup.(6), {testgroup.State, testgroup.Flicker}, 'model', 'interaction', 'varnames', {'State', 'Flicker'}); % data is column 6 - AMP 6/20/24
            disp('A Two-Way Anova was run because there are 2 IVs.');
            if  p(1,1) <= 0.05 % Anova for first IV
                t_results = sprintf( ' %s was significantly different across %s (F (%s,%s)= %s, p = %s). ', variablename, string(t{2,1}), string(t{2,3}), string(t{5,3}), string(t{2,6}), string(t{2,7}));
                disp(t_results);
            else
                t_results = sprintf('The results for %s were not significant. ', string(t{2,1})) ;
                disp(t_results);
            end

            if  p(2,1) <= 0.05 % Anova for second IV
                t_results = sprintf( ' %s was significantly different across %s (F (%s,%s)= %s, p = %s). ', variablename, string(t{3,1}), string(t{3,3}), string(t{5,3}), string(t{3,6}), string(t{3,7}));
                disp(t_results);
            else
                t_results = sprintf('The results for %s were not significant. ', string(t{3,1})) ;
                disp(t_results);
            end

            if  p(3,1) <= 0.05 % Anova for interaction
                t_results = sprintf( ' %s was significantly different across %s (F (%s,%s)= %s, p = %s). ', variablename, string(t{4,1}), string(t{4,3}), string(t{5,3}), string(t{4,6}), string(t{4,7}));
                disp(t_results);
            else
                t_results = sprintf('The results for %s were not significant. ', string(t{4,1})) ;
                disp(t_results);
            end

            manovafile = ([logpathname '/Stats/Anova_' variablename '_' testgroupname '_' date ]);
            save(manovafile,'t');
            % t=cell2table(t);
            text_multiflie = ([manovafile '.txt']);
            %writetable(t, text_multiflie,'Delimiter','\t','WriteRowNames',true);
            disp(t)

            % run post hoc comparisons between groups
            posthoc_comps = multcompare(stats, 'CriticalValueType', 'tukey-kramer'); %returns group, control group, lower limit, diff, upper limit, pval
            [results,~,~,gnames] = multcompare(stats, 'CType', 'tukey-kramer'); %returns posthoc paired comparisons test using correction for multiple comparisons
            posthocs_comps= array2table(results,"VariableNames", ["Group","Control Group", "Lower Limit", "Difference", "upper Limit", "P-value"]);
            posthocs_comps.("Group") = gnames(posthocs_comps.("Group")); % identify the group names
            posthocs_comps.("Control Group") = gnames(posthocs_comps.("Control Group")); % specify the comparsion groups column 2
            newmultifile = ([logpathname '/Stats/Anova_posthoc_comparisons_' variablename '_' region  '_' testgroupname '_' date]);
            save (newmultifile, 'posthocs_comps');
            text_multiflie = ([newmultifile '.txt']);
            writetable(posthocs_comps, text_multiflie,'Delimiter','\t','WriteRowNames',true);

        end
    end
    disp('the p value for the between groups test was' );
    % disp([ p ' for ' testgroupname 'scores on ' stringToBeFound ]);
    if p < .045
        disp ('woot! that is significant ');
      
    elseif p <= .05
        disp('it is just this side of significant..');
    else
        disp('it was not significant, sorry :(' );
    end
    disp('Between groups and post-hoc comparison stats were saved');
% end %prior end of sex iteration loop 

%% Create Graphs - a violin plot for each value
    disp(' ');
    disp('statistics and graphs will be displayed by cell unless edited');
    disp('making graphs now');

%% females and males graphed separately; females and males were exposed to different frequencies 6/13/24

    % "plot the Violin", somewhere around line 145, you can add 'Linestyle', 'none' to the "fill" function and that will remove the outline
    % figure;
    % updated_unblindedinfo.(3)=double(updated_unblindedinfo.(3));
    
    % set up here if different group orders based on sex or other variable
    % Males conditions in order as they appear on the graph
    male_graph_grouporder = {'CONTROL-NO STIMULATION', 'CONTROL-10HZ', 'STRESS-NO STIMULATION', 'STRESS-10HZ'};
    % Females conditions in order as they appear on the graph
    female_graph_grouporder = {'CONTROL-NO STIMULATION', 'CONTROL-40HZ', 'STRESS-NO STIMULATION', 'STRESS-40HZ'};

    % Customization for violin plots etc. below are optional for order, color,
    % edit for your specific experiment using decimal values by divide RGB/256
    color_ControlNoStim = [.45313 .45313 .45313];
    color_Control40Hz = [.81641 .81641 .81641];
    color_StressNoStim = [.52344 .10938 .19922];
    color_Stress40Hz = [.09766 .71484 .59766];
    color_Control10Hz = [.81641 .81641 .81641];
    color_Stress10Hz = [ .17578 .33594 .48438];

    g=cellstr(unique(testgroup.Condition)); % define unique experimental groups for the table

    %options ViolinColor, EdgeColor, MedianColor, BoxColor..
    figure('Position', [19 402 780 580]);    
    % figure('Position', [19 402 462 398]);
    
    if i == 1 % male data order
        testgroup_order = male_graph_grouporder;
    else % female data order
        testgroup_order = female_graph_grouporder;
    end
    % create the plot here, but it will edit the colors in the lines below
    vplot = violinplot(testgroup.(6), testgroup.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', testgroup_order);
    % this could probably be set up as an iteration through vplot ()
    % numbers and then adding 'color_' to the corresponding condition color
    if i == 1 % male data order and corresponding color order
        vplot(1).ViolinColor = color_ControlNoStim;
        vplot(2).ViolinColor = color_Control10Hz;
        vplot(3).ViolinColor = color_StressNoStim; 
        vplot(4).ViolinColor = color_Stress10Hz;

    else % female data order and corresponding colororder
        vplot(1).ViolinColor = color_ControlNoStim;
        vplot(2).ViolinColor = color_Control40Hz;
        vplot(3).ViolinColor = color_StressNoStim;
        vplot(4).ViolinColor = color_Stress40Hz;
    end
    
    % xticklabels('testgroup_order', 'Interpreter', 'none');
    xlabel('Condition');
    ylabel([variablename ' in ' region]);

    %Add significance bars using sigstar function
     %   * represents p<=0.05
    %  ** represents p<=1E-2
    % *** represents p<=1E-3

    for e = 1:height(posthocs_comps)
        pval = posthocs_comps{e,6};
        if pval<=0.05
            name1 = posthocs_comps{e,1};
            name1 = name1{1};

            % name1 = [beginning, name1];
            name1 = {name1};
            % name1 = name1(1);
            % name1 = name1(11:end);
            name2 = posthocs_comps{e,2};
            name2 = name2{1};
            % name2 = [beginning, name2];
            name2 = {name2};
            % name2 = name2(1};
            % name2 = name2(11:end);
            sigstar({[name1, name2]}, pval);
        end
    end
    saveas(gcf,[figpathname variablename '_' testgroupname '_violin_plot_' date '.png']);
    saveas(gcf,[figpathname variablename '_' testgroupname '_violin_plot_' date '.eps']);

    disp('')
    disp ('group graphs created and saved');
    
    %% If you want to graph the data by Filament ID that was averaged as a violin plot
    
    if avgcalc == 1
        % graph averages for subject and condition
        avgdata=grpstats(testgroup,["Condition","Subject","FilamentID"],["mean","sem"],"DataVars",variablename);
        % creage averages graph where each dot represents one animal
        figure('Position', [19 402 462 398]);
        vplot = violinplot(avgdata.(65), avgdata.Condition, 'BoxColor', [ 0 0 0], 'GroupOrder', {'WT No Flicker', '5XFAD No Flicker', '5XFAD 20Hz', '5XFAD 40Hz'});
        vplot(1).ViolinColor = color_WTNoFlicker;
        vplot(2).ViolinColor = color_5XFADNoFlicker;
        vplot(3).ViolinColor = color_5XFAD20Hz;
        vplot(4).ViolinColor = color_5XFAD40Hz;
        
        xlabel('Condition');
        ylabel(['Average ' variablename ' per cell in ' region]);
        saveas(gcf,[figpathname variablename '_cell_averages_violin_plot_' date '.png']);
    end
    disp('graphs have been saved');

    %% individual Violin Plots
    disp('making individual graphs now');
    %clf;
    %set up group order and color order
    % Males conditions in order as they appear on the graph (1)
    if i==1
    allConditions= {'CONTROL-NO STIMULATION', 'CONTROL-10HZ', 'STRESS-NO STIMULATION', 'STRESS-10HZ'};
    conditionColors = [color_ControlNoStim; color_Control10Hz; color_StressNoStim; color_Stress10Hz];

    % Females conditions in order as they appear on the graph (2)
    else
    allConditions = {'CONTROL-NO STIMULATION', 'CONTROL-40HZ', 'STRESS-NO STIMULATION', 'STRESS-40HZ'};
    conditionColors = [color_ControlNoStim; color_Control40Hz; color_StressNoStim; color_Stress40Hz];
    end

    % color per condition based on RGB values / 256 to convert to Matlab
    color_ControlNoStim = [.45313 .45313 .45313];
    color_Control40Hz = [.81641 .81641 .81641];
    color_StressNoStim = [.52344 .10938 .19922];
    color_Stress40Hz = [.09766 .71484 .59766];
    color_Control10Hz = [.81641 .81641 .81641];
    color_Stress10Hz = [.17578 .33594 .48438];

    % ## edited below here
    
    [X,Y] = ismember(testgroup.Condition,testgroup_order); % sort data based on testgroup order for graphing
    [~,Z] = sort(Y);
    testgroup=testgroup(Z,:); % sortconditions for testgroup
    plotSubjName = unique(cellstr(testgroup.Subject), 'stable'); % get unique sorted names keeping the order

%     [uniq_subj, ia] = unique(testgroup.Subject); % identify the unique subject name and index
%     uniq_cond = testgroup.Condition(ia); % identify matching conditions

    figure
    vplot = violinplot(testgroup.(6), testgroup.Subject, 'BoxColor', [ 0 0 0], 'GroupOrder', plotSubjName); %make violin plot
    hold on % pause to get an iteration

%%% edited above, below here is the same

    % sortrows(testgroup.Condition); % sort subject data based on condition
    % %sort rows based on condition
    plotSubjName = unique(cellstr(testgroup.Subject), 'stable'); %get sorted names keeping the order
    % figure % create figure area
    % vplot = violinplot(testgroup.(6), testgroup.Subject, 'BoxColor', [ 0 0 0], 'GroupOrder', plotSubjName); %make violin plot

    hold on
    [subj, subjI] = unique(testgroup.Subject, 'stable'); % subj = subject list
    subjCond = testgroup.Condition(subjI); % conditions for each subject
    subjColorMat = NaN(length(subj), 3); % create matrix of 3 column NaNs
    % 
    for c = 1:length(allConditions)
        isCond = strcmp(subjCond, allConditions{c});
        subjColorMat(isCond,:) = repmat(conditionColors(c,:), sum(isCond),1);
    end
    % % vplot = violinplot(testgroup.(6), testgroup.Subject, 'BoxColor', [ 0 0 0], 'GroupOrder', plotSubjName);
    % 
    for d = 1:size(subjColorMat,1)
        vplot(d).ViolinColor = [subjColorMat(d,:)];
        vplot(d).EdgeColor = [subjColorMat(d,:)];
    end
    
    xlabel('Subject ID')
    ylabel(variablename);
    hold off
    saveas(gcf,[figpathname variablename '_' testgroupname '_violin_plot_individual' date '.png']);
    saveas(gcf,[figpathname variablename '_' testgroupname '_violin_plot_individual' date '.eps']);


end 

disp('');
disp('*');
disp('*');
disp('*');
disp('Done.');