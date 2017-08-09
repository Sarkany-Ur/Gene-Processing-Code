%% Import data from spreadsheet

[~, ~, raw] = xlsread('Home Directory/Orthologs Between Species, Pared.xlsx','070517 Human Orthos');
raw = raw(2:end,1:11);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,4,5,8,11]);
raw = raw(:,[2,3,6,7,9,10]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));

% Create table
DRISC_Human_Orthologs = table;

% Allocate imported array to column variable names
DRISC_Human_Orthologs.MouseGeneName = cellVectors(:,1);
DRISC_Human_Orthologs.MouseGeneID = data(:,1);
DRISC_Human_Orthologs.MGIID = data(:,2);
% DRISC_Orthologs.MouseSymbol = cellVectors(:,2);
% DRISC_Orthologs.Species = cellVectors(:,3);
DRISC_Human_Orthologs.HumanGeneID = data(:,3);
DRISC_Human_Orthologs.HGNCID = data(:,4);
DRISC_Human_Orthologs.HumanSymbol = cellVectors(:,4);
DRISC_Human_Orthologs.DIOPTScore = data(:,5);
% DRISC_Orthologs.WeightedScore = data(:,6);
DRISC_Human_Orthologs.Rank = cellVectors(:,5);

% Clear temporary variables
clearvars data raw cellVectors R;

save('Home Directory/Disease Loaded Data/Data_DRISC_Human_Orthologs.mat')

clearvars DRISC_Human_Orthologs