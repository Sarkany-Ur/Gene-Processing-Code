%% Import data from spreadsheet
[~, ~, raw] = xlsread('Home Directory/Monarch - Gastrointestinal System Disease.xlsx','Monarch Gastrointestinal System');
raw = raw(:,[1:6,9:end]);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,5,6,7,8,9,10,11]);

% Create table
MONARCH_Disease_Gastrointestinal = table;

% Allocate imported array to column variable names
MONARCH_Disease_Gastrointestinal.NCBI_Gene_ID          = cellVectors(:,1);
MONARCH_Disease_Gastrointestinal.Gene_Name             = cellVectors(:,2);
% MONARCH_Abdomen.subject_taxon         = cellVectors(:,3);
% MONARCH_Abdomen.subject_taxon_label   = cellVectors(:,4);
MONARCH_Disease_Gastrointestinal.HPO_ID                = cellVectors(:,5);
MONARCH_Disease_Gastrointestinal.HPO_Phenotype         = cellVectors(:,6);
% MONARCH_Abdomen.evidence              = cellVectors(:,7);
% MONARCH_Abdomen.evidence_label        = cellVectors(:,8);
% MONARCH_Abdomen.source                = cellVectors(:,9);
% MONARCH_Abdomen.is_defined_by         = cellVectors(:,10);
% MONARCH_Abdomen.qualifier             = cellVectors(:,11);


%
% Clear temporary variables
clearvars raw cellVectors;

MONARCH_Disease_Gastrointestinal(1,:) = [];

save('Home Directory/Data_MONARCH_Disease_Gastrointestinal.mat')

clearvars MONARCH_Disease_Gastrointestinal