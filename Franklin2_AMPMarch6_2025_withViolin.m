%% instructions included above each section. run separately

%% A . loading and organizing data

% Run this section within the folder containing the excel files from
% imaris.

% This creates a structure called reconstructions with all parameters of
% interest for each cell. The data within this structure can be
% copied and pasted into other programs or made into it's own variable.
% It is primarily based on the file names, which are standardized across
% the dataset

reconstructions=[]
% addpath(genpath("C:\Users\cathe\Reconstructions\FlickerAstrocyteStudy")) %ADJUST TO lOCATION WHERE FILES ARE SAVED
addpath(genpath("T:\singer\Tina\Projects\Flicker Stress Mitigation_2025\IL")) %ADJUST TO lOCATION WHERE FILES ARE SAVED

%Open folder with desired files (i.e PL first, then IL)
%ACTIVE FOLDER DETERMINES WHICH REGION ANALYZED

listing=dir('**/*Filament.xls');  %list of every soma data file
files={listing.name};   %file names of every filament data file
   zooms=contains(files,'Zoom')
   files=files(~zooms) %remove zoomed images
start=1;    %to keep track of progress
reconstructions=[]

for i=1:length(files)
        filamentfile=files{i};
        spotsfile=strcat(filamentfile(1:end-12),['Spots.xls']);
        hullfile=strcat(filamentfile(1:end-12),['Hull.xls']);
        
        c=strsplit(filamentfile,'_');
        group=[c{1}];
        switch group;
            case 'Beta'
                sex='F';
                group='40H';
            case 'Delta'
                sex='F';
                group='CTL';
            case 'Gamma'
                sex='F';
                group='CUS';
            case 'Theta'
                sex='F';
                group='10H';
            case 'Alpha'
                sex='M';
                group='CTL';
            case 'Nu'
                sex='M';
                group='40H';
            case 'Pi'
                sex='M';
                group='CUS';
            case 'Zeta'
                sex='M';
                group='10H';
        end
        sample=[c{3}];
            sample=sample(1:2);
        side=[c{3}];
            side=side(end);
        location=[c{4} c{5}];
          
        soma=readmatrix(spotsfile,'Sheet','Overall');
        somanum=soma(1,2);
        hull=readmatrix(hullfile,'Sheet','Area');
        hullarea=hull(:,1);
        filament=readmatrix(filamentfile,'Sheet','Filament Dendrite Area (sum)');
        arbarea=filament(:,1);
        filament=readmatrix(filamentfile,'Sheet','Filament Dendrite Length (sum)');
        dendlength=filament(:,1);
        filament=readmatrix(filamentfile,'Sheet',28);
        junctionnum=filament(:,1); 
      
        reconstructions(start).group=group;
        reconstructions(start).sex=sex;
        reconstructions(start).location=location;
        reconstructions(start).sample=sample;
        reconstructions(start).side=side;
        reconstructions(start).somanum=somanum;
        reconstructions(start).hullarea=hullarea;
        reconstructions(start).arbarea=arbarea;
        reconstructions(start).dendlength=dendlength;
        reconstructions(start).junctionnum=junctionnum;
        
        start=start+1

        save('reconstructions','reconstructions');  %change filename as desired
        
end    


%% B . assign variables and create summary

% Section A can be run on all of the data (and the id of the groups are
% automatically associated with each set of measures) or can be run on a
% subset of the data. the following sections assume that Section A was run 
% on a subset of the data, distinguishing by cell type and brain region.
% i.e. astrocytes in the IL region. 

% Section B should be run for each region. After running, save the summary
% for later graphing using "summaryIL=summary"
% If you want the variables to include their source (i.e. instead of
% somanumCUS, you would have ILsomanumCUS), change name to 1 instead of 0
% and pre to the desire prefix

name = 1;
out = 1; 
pre='IL'; % CHANGE NAME TO MATCH REGION BEING ANALYZED

%define variables
somanumFCTL=[];somanumFCUS=[];somanumF10H=[];somanumF40H=[];somanumMCTL=[];somanumMCUS=[];somanumM10H=[];somanumM40H=[];
hullareaFCTL=[];hullareaFCUS=[];hullareaF10H=[];hullareaF40H=[];hullareaMCTL=[];hullareaMCUS=[];hullareaM10H=[];hullareaM40H=[];
arbareaFCTL=[];arbareaFCUS=[];arbareaF10H=[];arbareaF40H=[];arbareaMCTL=[];arbareaMCUS=[];arbareaM10H=[];arbareaM40H=[];
dendlengthFCTL=[];dendlengthFCUS=[];dendlengthF10H=[];dendlengthF40H=[];dendlengthMCTL=[];dendlengthMCUS=[];dendlengthM10H=[];dendlengthM40H=[];
junctionnumFCTL=[];junctionnumFCUS=[];junctionnumF10H=[];junctionnumF40H=[];junctionnumMCTL=[];junctionnumMCUS=[];junctionnumM10H=[];junctionnumM40H=[];
summary=[];

for i=1:length(reconstructions);
    switch [reconstructions(i).sex,reconstructions(i).group];
        case 'FCTL';
            somanumFCTL=[somanumFCTL;reconstructions(i).somanum];
            hullareaFCTL=[hullareaFCTL;reconstructions(i).hullarea];
            arbareaFCTL=[arbareaFCTL;reconstructions(i).arbarea];
            dendlengthFCTL=[dendlengthFCTL;reconstructions(i).dendlength];
            junctionnumFCTL=[junctionnumFCTL;reconstructions(i).junctionnum];
        case 'FCUS';
            somanumFCUS=[somanumFCUS;reconstructions(i).somanum];
            hullareaFCUS=[hullareaFCUS;reconstructions(i).hullarea];
            arbareaFCUS=[arbareaFCUS;reconstructions(i).arbarea];
            dendlengthFCUS=[dendlengthFCUS;reconstructions(i).dendlength];
            junctionnumFCUS=[junctionnumFCUS;reconstructions(i).junctionnum];
        case 'F10H';
            somanumF10H=[somanumF10H;reconstructions(i).somanum];
            hullareaF10H=[hullareaF10H;reconstructions(i).hullarea];
            arbareaF10H=[arbareaF10H;reconstructions(i).arbarea];
            dendlengthF10H=[dendlengthF10H;reconstructions(i).dendlength];
            junctionnumF10H=[junctionnumF10H;reconstructions(i).junctionnum];
        case 'F40H' ;
            somanumF40H=[somanumF40H;reconstructions(i).somanum];
            hullareaF40H=[hullareaF40H;reconstructions(i).hullarea];
            arbareaF40H=[arbareaF40H;reconstructions(i).arbarea];
            dendlengthF40H=[dendlengthF40H;reconstructions(i).dendlength];
            junctionnumF40H=[junctionnumF40H;reconstructions(i).junctionnum];
        case 'MCTL' ;
            somanumMCTL=[somanumMCTL;reconstructions(i).somanum];
            hullareaMCTL=[hullareaMCTL;reconstructions(i).hullarea];
            arbareaMCTL=[arbareaMCTL;reconstructions(i).arbarea];
            dendlengthMCTL=[dendlengthMCTL;reconstructions(i).dendlength];
            junctionnumMCTL=[junctionnumMCTL;reconstructions(i).junctionnum]; 
        case 'MCUS' ;
            somanumMCUS=[somanumMCUS;reconstructions(i).somanum];
            hullareaMCUS=[hullareaMCUS;reconstructions(i).hullarea];
            arbareaMCUS=[arbareaMCUS;reconstructions(i).arbarea];
            dendlengthMCUS=[dendlengthMCUS;reconstructions(i).dendlength];
            junctionnumMCUS=[junctionnumMCUS;reconstructions(i).junctionnum];  
        case 'M10H' ;
            somanumM10H=[somanumM10H;reconstructions(i).somanum];
            hullareaM10H=[hullareaM10H;reconstructions(i).hullarea];
            arbareaM10H=[arbareaM10H;reconstructions(i).arbarea];
            dendlengthM10H=[dendlengthM10H;reconstructions(i).dendlength];
            junctionnumM10H=[junctionnumM10H;reconstructions(i).junctionnum];
        case 'M40H' ;
            somanumM40H=[somanumM40H;reconstructions(i).somanum];
            hullareaM40H=[hullareaM40H;reconstructions(i).hullarea];
            arbareaM40H=[arbareaM40H;reconstructions(i).arbarea];
            dendlengthM40H=[dendlengthM40H;reconstructions(i).dendlength];
            junctionnumM40H=[junctionnumM40H;reconstructions(i).junctionnum];          
    end
end
summary=[mean(somanumFCTL) mean(somanumFCUS) mean(somanumF10H) mean(somanumF40H) mean(somanumMCTL) mean(somanumMCUS) mean(somanumM10H) mean(somanumM40H);
    std(somanumFCTL)/sqrt(length(somanumFCTL)) std(somanumFCUS)/sqrt(length(somanumFCUS)) std(somanumF10H)/sqrt(length(somanumF10H)) std(somanumF40H)/sqrt(length(somanumF40H)) std(somanumMCTL)/sqrt(length(somanumMCTL)) std(somanumMCUS)/sqrt(length(somanumMCUS)) std(somanumM10H)/sqrt(length(somanumM10H)) std(somanumM40H)/sqrt(length(somanumM40H));
    mean(arbareaFCTL) mean(arbareaFCUS) mean(arbareaF10H) mean(arbareaF40H) mean(arbareaMCTL) mean(arbareaMCUS) mean(arbareaM10H) mean(arbareaM40H);
    std(arbareaFCTL)/sqrt(length(arbareaFCTL)) std(arbareaFCUS)/sqrt(length(arbareaFCUS)) std(arbareaF10H)/sqrt(length(arbareaF10H)) std(arbareaF40H)/sqrt(length(arbareaF40H)) std(arbareaMCTL)/sqrt(length(arbareaMCTL)) std(arbareaMCUS)/sqrt(length(arbareaMCUS)) std(arbareaM10H)/sqrt(length(arbareaM10H)) std(arbareaM40H)/sqrt(length(arbareaM40H));
    mean(dendlengthFCTL) mean(dendlengthFCUS) mean(dendlengthF10H) mean(dendlengthF40H) mean(dendlengthMCTL) mean(dendlengthMCUS) mean(dendlengthM10H) mean(dendlengthM40H);
    std(dendlengthFCTL)/sqrt(length(dendlengthFCTL)) std(dendlengthFCUS)/sqrt(length(dendlengthFCUS)) std(dendlengthF10H)/sqrt(length(dendlengthF10H)) std(dendlengthF40H)/sqrt(length(dendlengthF40H)) std(dendlengthMCTL)/sqrt(length(dendlengthMCTL)) std(dendlengthMCUS)/sqrt(length(dendlengthMCUS)) std(dendlengthM10H)/sqrt(length(dendlengthM10H)) std(dendlengthM40H)/sqrt(length(dendlengthM40H));
    mean(junctionnumFCTL) mean(junctionnumFCUS) mean(junctionnumF10H) mean(junctionnumF40H) mean(junctionnumMCTL) mean(junctionnumMCUS) mean(junctionnumM10H) mean(junctionnumM40H);
    std(junctionnumFCTL)/sqrt(length(junctionnumFCTL)) std(junctionnumFCUS)/sqrt(length(junctionnumFCUS)) std(junctionnumF10H)/sqrt(length(junctionnumF10H)) std(junctionnumF40H)/sqrt(length(junctionnumF40H)) std(junctionnumMCTL)/sqrt(length(junctionnumMCTL)) std(junctionnumMCUS)/sqrt(length(junctionnumMCUS)) std(junctionnumM10H)/sqrt(length(junctionnumM10H)) std(junctionnumM40H)/sqrt(length(junctionnumM40H));
    mean(hullareaFCTL) mean(hullareaFCUS) mean(hullareaF10H) mean(hullareaF40H) mean(hullareaMCTL) mean(hullareaMCUS) mean(hullareaM10H) mean(hullareaM40H);
    std(hullareaFCTL)/sqrt(length(hullareaFCTL)) std(hullareaFCUS)/sqrt(length(hullareaFCUS)) std(hullareaF10H)/sqrt(length(hullareaF10H)) std(hullareaF40H)/sqrt(length(hullareaF40H)) std(hullareaMCTL)/sqrt(length(hullareaMCTL)) std(hullareaMCUS)/sqrt(length(hullareaMCUS)) std(hullareaM10H)/sqrt(length(hullareaM10H)) std(hullareaM40H)/sqrt(length(hullareaM40H))];

% outliers, if desired:
if out
    base_names = {'somanumFCTL', 'somanumFCUS', 'somanumF10H', 'somanumF40H', 'somanumMCTL', 'somanumMCUS', 'somanumM10H', 'somanumM40H', 'hullareaFCTL', 'hullareaFCUS', 'hullareaF10H', 'hullareaF40H', 'hullareaMCTL', 'hullareaMCUS', 'hullareaM10H', 'hullareaM40H', 'arbareaFCTL', 'arbareaFCUS', 'arbareaF10H', 'arbareaF40H', 'arbareaMCTL', 'arbareaMCUS', 'arbareaM10H', 'arbareaM40H', 'dendlengthFCTL', 'dendlengthFCUS', 'dendlengthF10H', 'dendlengthF40H', 'dendlengthMCTL', 'dendlengthMCUS', 'dendlengthM10H', 'dendlengthM40H', 'junctionnumFCTL', 'junctionnumFCUS', 'junctionnumF10H', 'junctionnumF40H', 'junctionnumMCTL', 'junctionnumMCUS', 'junctionnumM10H', 'junctionnumM40H'};
    for i=1:length(base_names)
        var_name=base_names{i};
        data = evalin('base', var_name);
        outliers = isoutlier(data, 'grubbs');
        cleaned_data = data(~outliers);
        assignin('base', var_name, cleaned_data);
    end
end

% renaming, if desired: 
if name
    base_names = {'summary','somanumFCTL', 'somanumFCUS', 'somanumF10H', 'somanumF40H', 'somanumMCTL', 'somanumMCUS', 'somanumM10H', 'somanumM40H', 'hullareaFCTL', 'hullareaFCUS', 'hullareaF10H', 'hullareaF40H', 'hullareaMCTL', 'hullareaMCUS', 'hullareaM10H', 'hullareaM40H', 'arbareaFCTL', 'arbareaFCUS', 'arbareaF10H', 'arbareaF40H', 'arbareaMCTL', 'arbareaMCUS', 'arbareaM10H', 'arbareaM40H', 'dendlengthFCTL', 'dendlengthFCUS', 'dendlengthF10H', 'dendlengthF40H', 'dendlengthMCTL', 'dendlengthMCUS', 'dendlengthM10H', 'dendlengthM40H', 'junctionnumFCTL', 'junctionnumFCUS', 'junctionnumF10H', 'junctionnumF40H', 'junctionnumMCTL', 'junctionnumMCUS', 'junctionnumM10H', 'junctionnumM40H'};
    for i = 1:length(base_names)
    new_var_name = [base_names{i},pre];
    eval([new_var_name, ' = [];']);
    eval([new_var_name, ' = ', base_names{i}, ';']);
    end
end

%% Graphing - bar charts
% now you have all the variables in your workspace. use the following code
% to create bar charts by inputing the names of the variables you wish to
% plot into this function. 

function pltCAT(vars,y_limit)
    % vars: Cell array of variable names to be plotted (e.g., {'somanumCUSPLA','somanumCUSILA'})
    % y_limit (optional): A two-element vector specifying the y-axis range, e.g., [0, 100]. 
    
% Under the "%Customize plot" line, make adjustments to visuals and labels.
    %Customize plot
    plot_size=[3,3]
    font_size=13
    error_bar_width=3
    title_text={'Male', 'Arborization Area'} % CHANGE TO MATCH SEX
    x_tick_labels={'Cntl','NoStim','10Hz','40Hz'}
    y_label='Arborization Area (\mum^2)' % CHANGE TO MATCH STATISTIC
    bar_colors = {[116/255 116/255 116/255],[134/255 28/255 51/255],[45/255 86/255 124/255],[25/255 183/255 153/255]}
    error_bar_color = {[116/255 116/255 116/255],[134/255 28/255 51/255],[45/255 86/255 124/255],[25/255 183/255 153/255]}

    % Initialize arrays to hold means and standard errors
    means = zeros(1,length(vars));
    sems = zeros(1,length(vars));
    data_all = cell(1, length(vars));
    % Loop over the variables and calculate the mean and standard error
    for i = 1:length(vars)
        data=evalin('base',vars{i});
        means(i)=mean(data,'omitnan');
        sems(i)=std(data,'omitnan')/sqrt(length(data));
        data_all{i}=data;
    end

    % Create the bar plot
    figure('Position', [100, 100, plot_size(1)*100, plot_size(2)*100]); % Set the figure size
    % b=bar(means, 'FaceColor', 'flat', 'EdgeColor','none'); 
    x_positions = linspace(1, length(vars), length(vars)); % Evenly distribute bars
    b = bar(x_positions, means, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.8); % Keep normal width
    xlim([min(x_positions) - 0.5, max(x_positions) + 0.5]);
    hold on;
    for i = 1:length(vars)
        b.CData(i, :) = bar_colors{i};
    end

    % % % Add error bars, if not including the data points themselves (chunk
    % of code below)
    % for i = 1:length(means)
    %     errorbar(i, means(i), sems(i), 'Color', error_bar_color{i}, 'LineWidth',error_bar_width,'LineStyle', 'none'); 
    % end

    for i = 1:length(vars)
        if ~isempty(data_all{i})
            data = data_all{i};
            % Compute density estimation using ksdensity
            [density, y_values] = ksdensity(data); % Kernel density estimation
            % Normalize density to scale jitter width
            density = density / max(density); % Scale between 0 and 1
            max_jitter = 0.18; % Adjust this for overall width of jitter
            % Map data points to corresponding density-based jitter
            x_jitter = zeros(size(data));
            for j = 1:length(data)
                % Find closest y-value in density estimate
                [~, idx] = min(abs(y_values - data(j)));
                % Scale jitter width based on density
                x_jitter(j) = (rand - 0.5) * density(idx) * max_jitter * 2; 
            end
            % Plot jittered scatter points
            scatter(i + x_jitter, data, 25, bar_colors{i}, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerFaceAlpha', 0.6);
        end
    end

    if nargin > 1 && ~isempty(y_limit)
        ylim(y_limit);  % Set the y-axis limits manually if y_limit is provided
    end

    % Aesthetics
    title(title_text, 'FontSize', font_size);
    ylabel(y_label, 'FontSize', font_size);
    xticks(x_positions);
    xticklabels(x_tick_labels);
    box off;
    set(gca, 'Color', 'none');
    set(gca,'tickdir','out');
    set(gca, 'FontSize', font_size);
    set(gca, 'FontName', 'Calibri')
    set(gca, 'XTickLabelRotation', 0);  % Rotate x-tick labels if needed
    hold off;

end

%run with a command like this line - alter here to pick your groups
pltCAT({'arbareaMCTL', 'arbareaMCUS', 'arbareaM10H','arbareaM40H'},[0 2500])


%% %% G. Violin plots
% run separately for PL and IL (remember to change graph title)
% generates violin plots for both male and female animals

reconstructionsM=[];
reconstructionsF=[];
for i = 1:length(reconstructions)
    if strcmp(reconstructions(i).sex, 'M')  % Check if sex is 'M' for male
        reconstructionsM = [reconstructionsM, reconstructions(i)];  % Add to male structure
    elseif strcmp(reconstructions(i).sex, 'F')  % Check if sex is 'F' for female
        reconstructionsF = [reconstructionsF, reconstructions(i)];  % Add to female structure
    end
end

%Male
groupColors = struct('CTL', [116/255 116/255 116/255], ...
                     'CUS', [134/255 28/255 51/255], ...
                     'G10H', [45/255 86/255 124/255], ...
                     'G40H', [25/255 183/255 153/255]);
samples=unique({reconstructionsM.sample});
groups=unique({reconstructionsM.group});
groups = regexprep(groups, '^(10H|40H)$', 'G$1');
data = cell(length(samples), 1);
groupLabels = cell(length(samples), 1);
groupColorsMapped = zeros(length(samples), 3);
for i=1:length(samples)
    sampleIDx=strcmp({reconstructionsM.sample}, samples{i});
    data{i} = vertcat(reconstructionsM(sampleIDx).arbarea);
    groupName = reconstructionsM(find(sampleIDx, 1)).group;
        if any(strcmp(groupName, {'10H', '40H'}))
            groupName = ['G' groupName];
        end
    groupLabels{i} = groupName; 
    groupColorsMapped(i, :) = groupColors.(groupLabels{i});
end
groupOrder = {'CTL', 'CUS', 'G10H', 'G40H'};
[~, sortIdx] = ismember(groupLabels, groupOrder);
[~, reorderIdx] = sort(sortIdx);
groupLabels = groupLabels(reorderIdx); 
data = data(reorderIdx);
groupColorsMapped = groupColorsMapped(reorderIdx, :);
samples=samples(reorderIdx);
for i=1:length(data)
    outliers=isoutlier(data{i},'grubbs');
    data{i}=data{i}(~outliers);
end
figure;
hold on;
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 5, 3]);
for i = 1:length(samples)
    [f, xi] = ksdensity(data{i});
    fill([-f + i; f + i], [xi; flipud(xi)], groupColorsMapped(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    density = f / max(f); 
    max_jitter = 0.2;
    x_jitter = zeros(size(data{i}));
    for j = 1:length(data{i})
        [~, idx] = min(abs(xi - data{i}(j)));
        x_jitter(j) = (rand - 0.5) * density(idx) * max_jitter * 2; 
    end
        scatter(i + x_jitter, data{i}, 25, groupColorsMapped(i,:), 'filled','MarkerFaceAlpha', 0.6);
end
xlim([0 length(samples)+1]);
ylim([0,inf]);
ylabel('Arborization Area (\mum^2)');
set(gca,'XTickLabel',[]);
set(gca, 'XTick', [])
title('IL Astrocyte Arborization Area by Animal (Male)'); % change name to match region
set(gca, 'Color', 'none');
set(gca,'tickdir','out');
set(gca, 'FontSize', 12);
set(gca, 'FontName', 'Calibri')
legendLabels={'Cntl','NoStim','10Hz','40Hz'};
legendHandles = zeros(1, length(groupOrder));
for i = 1:length(groupOrder)
    legendHandles(i) = patch(nan, nan, 'w', 'FaceColor', groupColors.(groupOrder{i}), 'EdgeColor', 'none');
end
lgd=legend(legendHandles, legendLabels);
lgd.ItemTokenSize = [8, 8];
lgd.Box='off';
lgd.FontSize=9
box off;
hold off;

%repeated for females
groupColors = struct('CTL', [116/255 116/255 116/255], ...
                     'CUS', [134/255 28/255 51/255], ...
                     'G10H', [45/255 86/255 124/255], ...
                     'G40H', [25/255 183/255 153/255]);
samples=unique({reconstructionsF.sample});
groups=unique({reconstructionsF.group});
groups = regexprep(groups, '^(10H|40H)$', 'G$1');
data = cell(length(samples), 1);
groupLabels = cell(length(samples), 1);
groupColorsMapped = zeros(length(samples), 3);
for i=1:length(samples)
    sampleIDx=strcmp({reconstructionsF.sample}, samples{i});
    data{i} = vertcat(reconstructionsF(sampleIDx).arbarea);
    groupName = reconstructionsF(find(sampleIDx, 1)).group;
        if any(strcmp(groupName, {'10H', '40H'}))
            groupName = ['G' groupName];
        end
    groupLabels{i} = groupName; 
    groupColorsMapped(i, :) = groupColors.(groupLabels{i});
end

groupOrder = {'CTL', 'CUS', 'G10H', 'G40H'};
[~, sortIdx] = ismember(groupLabels, groupOrder);
[~, reorderIdx] = sort(sortIdx);
groupLabels = groupLabels(reorderIdx); 
data = data(reorderIdx);
groupColorsMapped = groupColorsMapped(reorderIdx, :);
samples=samples(reorderIdx);
for i=1:length(data)
    outliers=isoutlier(data{i},'grubbs');
    data{i}=data{i}(~outliers);
end


figure;
hold on;
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 5, 3]);
for i = 1:length(samples)
    [f, xi] = ksdensity(data{i});
    fill([-f + i; f + i], [xi; flipud(xi)], groupColorsMapped(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    density = f / max(f); 
    max_jitter = 0.2;
    x_jitter = zeros(size(data{i}));
    for j = 1:length(data{i})
        [~, idx] = min(abs(xi - data{i}(j)));
        x_jitter(j) = (rand - 0.5) * density(idx) * max_jitter * 2; 
    end
        scatter(i + x_jitter, data{i}, 25, groupColorsMapped(i,:), 'filled','MarkerFaceAlpha', 0.6);
end
xlim([0 length(samples)+1]);
ylim([0,inf]);
ylabel('Arborization Area (\mum^2)');
set(gca,'XTickLabel',[]);
set(gca, 'XTick', [])
title('IL Astrocyte Arborization Area by Animal (Female)'); % change name to match region
set(gca, 'Color', 'none');
set(gca,'tickdir','out');
set(gca, 'FontSize', 12);
set(gca, 'FontName', 'Calibri')
legendLabels={'Cntl','NoStim','10Hz','40Hz'};
legendHandles = zeros(1, length(groupOrder));
for i = 1:length(groupOrder)
    legendHandles(i) = patch(nan, nan, 'w', 'FaceColor', groupColors.(groupOrder{i}), 'EdgeColor', 'none');
end
lgd=legend(legendHandles, legendLabels);
lgd.ItemTokenSize = [8, 8];
lgd.Box='off';
lgd.FontSize=9
box off;
hold off;

%% Stats - 1 way anova and post hoc

function runStatsTests(vars)
    % vars: Cell array of variable names to be tested (e.g., {'dendlengthFCTL', 'dendlengthFCUS', 'dendlengthF10H','dendlengthF40H'})
    % This function will run ANOVA, t-tests, and Tukey's post-hoc tests.

    % Check if the variables exist in the workspace
    data = cell(1, length(vars));
    for i = 1:length(vars)
        if evalin('base', ['exist(''', vars{i}, ''', ''var'')'])
            % Get the current variable data
            data{i} = evalin('base', vars{i});
        else
            warning(['Variable ', vars{i}, ' does not exist in the workspace.']);
            data{i} = NaN; % Assign NaN if the variable does not exist
        end
    end

    data_matrix = cell2mat(data(:));
    group_labels = [];
    for i = 1:length(vars)
        group_labels = [group_labels; repmat(i, length(data{i}), 1)];
    end

    % Perform ANOVA (one-way)
    [p_anova, tbl, stats] = anova1(data_matrix,group_labels);

    % Display ANOVA result
    fprintf('ANOVA p-value: %.4f\n', p_anova);

    % If ANOVA is significant, run post hoc test (Tukey's HSD)
    if p_anova < 0.05
        [c, m, h, gnames] = multcompare(stats);
        disp(c);  % Display post hoc results (group comparisons)
    else
        fprintf('ANOVA is not significant. No post-hoc test is performed.\n');
    end

    % Run t-tests if there are exactly two variables
    if length(vars) == 2
        fprintf('\nRunning t-test between %s and %s...\n', vars{1}, vars{2});
        [h_ttest, p_ttest] = ttest2(data{1}, data{2});
        fprintf('t-test p-value: %.4f\n', p_ttest);
    end
end

% run this line to use the function (adjust inputs as necessary)
runStatsTests({'arbareaFCTLIL', 'arbareaFCUSIL', 'arbareaF10HIL','arbareaF40HIL'})


%% Stats - 2-way anova and post hoc

function run2FactorANOVA(vars, factor1_labels, factor2_labels)
    % vars: Cell array of variable names to be plotted (e.g., {'dendlengthFCUS', 'dendlengthF10H','dendlengthF40H','dendlengthMCUS', 'dendlengthM10H','dendlengthM40H'})
    % factor1_labels: Cell array of group labels for the first factor, like sex (e.g., {'F', 'F', 'F','M','M','M'})
    % factor2_labels: Cell array of group labels for the second factor, like flicker (e.g., {'00H', '10H', '40H', '00H','10H','40H'})

    % Initialize a cell array to hold the data for each group
    data = cell(1, length(vars));

    % Loop through the variable names to gather the data from the workspace
    for i = 1:length(vars)
        % Check if the variable exists in the workspace
        if evalin('base', ['exist(''', vars{i}, ''', ''var'')'])
            % Get the current variable data
            data{i} = evalin('base', vars{i});
        else
            warning(['Variable ', vars{i}, ' does not exist in the workspace.']);
            data{i} = NaN; % Assign NaN if the variable doesn't exist
        end
    end

    % Ensure that all variables are column vectors (XXX by 1)
    for i = 1:length(vars)
        if size(data{i}, 2) ~= 1
            error(['The variable ', vars{i}, ' is not a column vector.']);
        end
    end

    % Check if the number of factor labels matches the number of variables
    if length(factor1_labels) ~= length(vars) || length(factor2_labels) ~= length(vars)
        error('The number of factor labels does not match the number of variables.');
    end

    % Concatenate the data from each variable into a single matrix (data_matrix)
    data_matrix = cell2mat(data(:));

    % Create a group vector indicating which group each data point belongs to
    group_labels_1 = [];
    group_labels_2 = [];

    % Create group labels for each factor
    for i = 1:length(vars)
        % Factor 1: Group ID (e.g., CTL, CUS)
        group_labels_1 = [group_labels_1; repmat(factor1_labels{i}, length(data{i}), 1)];
        % Factor 2: Condition ID (e.g., IL, PL)
        group_labels_2 = [group_labels_2; repmat(factor2_labels{i}, length(data{i}), 1)];
    end

    % Combine the two factors into a single factor (interaction term)
    % Create an interaction term from both factors
    interaction_labels = strcat(num2str(group_labels_1), '-', num2str(group_labels_2));

    % Perform multi-factor ANOVA with interaction
    % Combine all factors into a single factor label array
    all_labels = [group_labels_1, group_labels_2];

    % Create a table for the ANOVA
    tbl = table(data_matrix, group_labels_1, group_labels_2, 'VariableNames', {'Data', 'Factor1', 'Factor2'});

    % Perform two-way ANOVA
    [p_anova, tbl, stats] = anovan(tbl.Data, {tbl.Factor1, tbl.Factor2}, 'model', 'interaction', 'varnames', {'Factor1', 'Factor2'});

    % Display the ANOVA table
    disp('Two-Way ANOVA Results:');
    disp(tbl);

    % Perform post-hoc test (Tukey's HSD) to compare means between different conditions
    disp('Post-hoc Tukey test results:');
    %multcompare(stats, 'Dimension', [1 2]);
    comparison_results = multcompare(stats, 'Dimension', [1 2]);
    disp('Pairwise Comparison p-values:');
    for i = 1:size(comparison_results, 1)
        fprintf('Group %d vs Group %d: p = %.4f\n', ...
            comparison_results(i, 1), comparison_results(i, 2), comparison_results(i, 6));
    end
end

% run this line to use function, changing groups as desired
runMultiFactorANOVA({'dendlengthFCUS', 'dendlengthF10H','dendlengthF40H','dendlengthMCUS', 'dendlengthM10H','dendlengthM40H'}, {'F', 'F', 'F','M','M','M'}, {'00H', '10H', '40H', '00H','10H','40H'})

%% stats - multifactor anova and stats

function runMultiFactorANOVA(vars, factor1_labels, factor2_labels, varargin)
    % vars: Cell array of variable names to be analyzed
    % factor1_labels: Cell array of group labels for the first factor
    % factor2_labels: Cell array of group labels for the second factor
    % varargin: Additional factors (e.g., factor3_labels, factor4_labels, ...)

    numFactors = 2 + length(varargin);  % Determine the number of factors
    factorLabels = [{factor1_labels}, {factor2_labels}, varargin];  % Store all factors

    % Initialize cell array for data storage
    data = cell(1, length(vars));

    % Loop through the variable names to gather the data from the workspace
    for i = 1:length(vars)
        if evalin('base', ['exist(''', vars{i}, ''', ''var'')'])
            data{i} = evalin('base', vars{i});
        else
            warning(['Variable ', vars{i}, ' does not exist in the workspace.']);
            data{i} = NaN;
        end
    end

    % Ensure all variables are column vectors
    for i = 1:length(vars)
        if size(data{i}, 2) ~= 1
            error(['The variable ', vars{i}, ' is not a column vector.']);
        end
    end

    % Check if the number of factor labels matches the number of variables
    for f = 1:numFactors
        if length(factorLabels{f}) ~= length(vars)
            error(['The number of labels for Factor ', num2str(f), ' does not match the number of variables.']);
        end
    end

    % Concatenate the data into a single matrix
    data_matrix = cell2mat(data(:));

    % Create group vectors for each factor
    group_labels = cell(numFactors, 1);
    for f = 1:numFactors
        group_labels{f} = {}; % Initialize as a cell array
    end

    for i = 1:length(vars)
        for f = 1:numFactors
            group_labels{f} = [group_labels{f}; repmat({factorLabels{f}{i}}, length(data{i}), 1)];
        end
    end


    % Convert all factors to categorical for ANOVA
    for f = 1:numFactors
        group_labels{f} = categorical(group_labels{f});
    end

    % Create a table for the ANOVA
    factorNames = strcat("Factor", string(1:numFactors));
    tbl = table(data_matrix, group_labels{:}, 'VariableNames', ['Data', factorNames]);

    % Perform multi-factor ANOVA
    modelSpec = 'interaction';  % Includes all interactions
    [p_anova, anova_table, stats] = anovan(data_matrix, group_labels, 'model', modelSpec, 'varnames', factorNames);


    % Display the ANOVA table
    disp(['Factorial ANOVA (', num2str(numFactors), '-way) Results:']);
    disp(anova_table);

    % Perform post-hoc Tukey test for all factors
    disp('Post-hoc Tukey test results:');
    for f = 1:numFactors
        disp(['Post-hoc comparisons for ', factorNames(f), ':']);
        comparison_results = multcompare(stats, 'Dimension', f);
        for i = 1:size(comparison_results, 1)
            fprintf('Group %d vs Group %d: p = %.4f\n', ...
                comparison_results(i, 1), comparison_results(i, 2), comparison_results(i, 6));
        end
    end
end

%% AMP added on March 6 to extract data and save structures in .csv
updatedfile=strcat(pwd, '/Data/', location, '_morphology_data__appended_', date, '.mat'); % saves out reconstruction matlab structure of doubles
save(updatedfile,'reconstructions'); % save the updated data structure under Data folder
% just for soma data below
all_data = struct2table(reconstructions); % convert reconstructions to a matlab table
soma_data = all_data(:,1:6); % select soma range
updatedfile=strcat(pwd, '/Data/', location, '_somanum_morphology_data__appended_', date, '.mat'); % saves out reconstruction matlab structure
save(updatedfile,'soma_data'); % save the updated data structure under Data folder
updated_csvfile=strcat(pwd, '/Data/', location, '_somanum_morphology_data_appended_', date, '.csv');
writetable(soma_data, updated_csvfile);

for j = 7:10 % for the remaining morphology columns 7-10
    aggregate = table; % create blank table
    var = string(all_data.Properties.VariableNames(j)); % find variable name
    for k = 1:height(all_data)
        rowval = (all_data(k,1:5)); % iterate through each row
        sup_data = table(all_data.(j){k,1}); % find the data for that sample
        repval = height(sup_data); % how much data per sample
        dummyvals = repmat(rowval,repval,1); % create a table of identifying info
        newdata = [dummyvals sup_data]; % combine the identifiers with variable data for that sample
        newdata.Properties.VariableNames(6) = var; % relable columns
        aggregate = [aggregate; newdata]; % combine the data with the next animal's data
    end
    updatedfile=strcat(pwd, '/Data/', location, '_', var, '_morphology_data__appended_', date, '.mat'); % saves out reconstruction matlab structure of doubles
    save(updatedfile,'aggregate'); % save the updated data structure under Data folder
    updated_csvfile=strcat(pwd, '/Data/', location, '_', var, '_morphology_data_', '_appended_', date, '.csv');
    writetable(aggregate, updated_csvfile);
end
disp ('The Unblinded data was saved as matlab structure and csv file under the Data path for the Parameters chosen.');

