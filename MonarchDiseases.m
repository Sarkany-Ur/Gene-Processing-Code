
%% ------------------------ Processing Script ------------------------
%% Load the datasets
% July 1, 2017

% Set up the directory paths and load data
clear all; clc
CurrentDirectory = pwd;
Data_File_Folder  = 'Home Directory';
cd(Data_File_Folder)

% The following is for finding a human gene and its related diseases
% Load a human-relevant data set
load('Data_MONARCH_Disease_Cardiovascular.mat')
load('Data_MONARCH_Disease_Endocrine.mat')
load('Data_MONARCH_Disease_Gastrointestinal.mat')
load('Data_MONARCH_Disease_Immune.mat')
load('Data_MONARCH_Disease_Integument.mat')
load('Data_MONARCH_Disease_Musculoskeletal.mat')
load('Data_MONARCH_Disease_Nervous_System.mat')
load('Data_MONARCH_Disease_Reproductive_System.mat')
load('Data_MONARCH_Disease_Respiratory_System.mat')
load('Data_MONARCH_Disease_Thoracic.mat')
load('Data_MONARCH_Disease_Urinary.mat')
load('Mouse_Gene_List.mat')
load('Data_DRISC_Human_Orthologs.mat');

cd(CurrentDirectory)
clear Data_File_Folder CurrentDirectory

% -------------------------------------------------------------------------
% DRISC Mouse-Human Ortho Data
% -------------------------------------------------------------------------
DRISC_Human_Orthologs       = table2cell(DRISC_Human_Orthologs);
Array_Length_DRISC_Human    = max(size(DRISC_Human_Orthologs));
Pared_DRISC_Human_Orthologs(Array_Length_DRISC_Human) = struct();

for i = 1:Array_Length_DRISC_Human
    Pared_DRISC_Human_Orthologs(i).Mouse_Gene_Name    = char(DRISC_Human_Orthologs(i,1));
    Pared_DRISC_Human_Orthologs(i).Mouse_Gene_ID      = cell2mat(DRISC_Human_Orthologs(i,2));
    Pared_DRISC_Human_Orthologs(i).MGI_ID             = cell2mat(DRISC_Human_Orthologs(i,3));
    Pared_DRISC_Human_Orthologs(i).Human_Gene_ID      = cell2mat(DRISC_Human_Orthologs(i,4));
    Pared_DRISC_Human_Orthologs(i).HGN_ID             = cell2mat(DRISC_Human_Orthologs(i,5));
    Pared_DRISC_Human_Orthologs(i).Human_Gene_Name    = char(DRISC_Human_Orthologs(i,6));
    Pared_DRISC_Human_Orthologs(i).DIOPT_Score        = cell2mat(DRISC_Human_Orthologs(i,7));
    Pared_DRISC_Human_Orthologs(i).Rank               = char(DRISC_Human_Orthologs(i,8));
end

% Run each data set through a small paring function to parse out parts we are interested in.
Pared_MONARCH_Disease_Cardiovascular = Data_Paring(MONARCH_Disease_Cardiovascular);
Array_Length_Cardiovascular = size(Pared_MONARCH_Disease_Cardiovascular,2);           % Find the length of the data set

Pared_MONARCH_Disease_Endocrine = Data_Paring(MONARCH_Disease_Endocrine);
Array_Length_Endocrine = size(Pared_MONARCH_Disease_Endocrine,2);           % Find the length of the data set

Pared_MONARCH_Disease_Gastrointestinal = Data_Paring(MONARCH_Disease_Gastrointestinal);
Array_Length_Gastrointestinal = size(Pared_MONARCH_Disease_Gastrointestinal,2);           % Find the length of the data set

Pared_MONARCH_Disease_Immune = Data_Paring(MONARCH_Disease_Immune);
Array_Length_Immune = size(Pared_MONARCH_Disease_Immune,2);           % Find the length of the data set

Pared_MONARCH_Disease_Integument = Data_Paring(MONARCH_Disease_Integument);
Array_Length_Integument = size(Pared_MONARCH_Disease_Integument,2);           % Find the length of the data set

Pared_MONARCH_Disease_Musculoskeletal = Data_Paring(MONARCH_Disease_Musculoskeletal);
Array_Length_Musculoskeletal = size(Pared_MONARCH_Disease_Musculoskeletal,2);           % Find the length of the data set

Pared_MONARCH_Disease_Nervous_System = Data_Paring(MONARCH_Disease_Nervous_System);
Array_Length_Nervous_System = size(Pared_MONARCH_Disease_Nervous_System,2);           % Find the length of the data set

Pared_MONARCH_Disease_Reproductive_System = Data_Paring(MONARCH_Disease_Reproductive_System);
Array_Length_Reproductive_System = size(Pared_MONARCH_Disease_Reproductive_System,2);           % Find the length of the data set

Pared_MONARCH_Disease_Respiratory_System = Data_Paring(MONARCH_Disease_Respiratory_System);
Array_Length_Respiratory_System = size(Pared_MONARCH_Disease_Respiratory_System,2);           % Find the length of the data set

Pared_MONARCH_Disease_Thoracic = Data_Paring(MONARCH_Disease_Thoracic);
Array_Length_Thoracic = size(Pared_MONARCH_Disease_Thoracic,2);           % Find the length of the data set

Pared_MONARCH_Disease_Urinary = Data_Paring(MONARCH_Disease_Urinary);
Array_Length_Urinary = size(Pared_MONARCH_Disease_Urinary,2);           % Find the length of the data set


%% Disease Search - Use for finding Human diseases

for x = 1:size(Mouse_Gene_List,1)
    
    Search_Term = Mouse_Gene_List{x,1};
    
    % -------------------------------------------------------------------------
    % Human Gene Search
    % -------------------------------------------------------------------------
    Indices_Human    = zeros(Array_Length_DRISC_Human,1);
    for i = 1:Array_Length_DRISC_Human
        Indices_Human(i) = strcmp(Pared_DRISC_Human_Orthologs(i).Mouse_Gene_Name, Search_Term);
    end
    Indices_Human = logical(Indices_Human);
    Disease_Results_Human = struct2cell(Pared_DRISC_Human_Orthologs(Indices_Human));
    
    Gene_Comparison_Human = cell(size(Disease_Results_Human,3),23);
    
    for i = 1:(size(Disease_Results_Human,3))
        Gene_Comparison_Human(i,1) = squeeze(Disease_Results_Human(1,:,i));     % Mouse Gene
        Gene_Comparison_Human(i,2) = squeeze(Disease_Results_Human(6,:,i));     % Human Gene
    end
    
    Gene_Comparison_Human(1,1)
    size(Disease_Results_Human,3);
    sort(Gene_Comparison_Human(:,2))';
    
    
    % Load specialized data sets and process them
    % Search the disease database by gene name, find the number of associated diseases
    % -------------------------------------------------------------------------
    Store_Index = 3;
     
    % Look at this first one for comments on additional useful output
    for i = 1:size(Gene_Comparison_Human,1)
        
        % Remove ";" to see the term it is currently looking at
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Cardiovascular,1);
        for j = 1:Array_Length_Cardiovascular
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Cardiovascular(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Cardiovascular(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        % Remove ";" from these to see the ACTUAL PHENOTYPES for each searched-for term
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Endocrine,1);
        
        for j = 1:Array_Length_Endocrine
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Endocrine(j).Gene_Name, Search_Term); % Find where the term is
        end
        
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Endocrine(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Gastrointestinal,1);
        
        for j = 1:Array_Length_Gastrointestinal
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Gastrointestinal(j).Gene_Name, Search_Term); % Find where the term is
        end
        
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Gastrointestinal(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Immune,1);
        for j = 1:Array_Length_Immune
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Immune(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Immune(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Integument,1);
        for j = 1:Array_Length_Integument
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Integument(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Integument(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Musculoskeletal,1);
        for j = 1:Array_Length_Musculoskeletal
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Musculoskeletal(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Musculoskeletal(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Nervous_System,1);
        for j = 1:Array_Length_Nervous_System
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Nervous_System(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Nervous_System(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human)
        size(Disease_Results_Human,3)
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Reproductive_System,1);
        for j = 1:Array_Length_Reproductive_System
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Reproductive_System(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Reproductive_System(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Respiratory_System,1);
        for j = 1:Array_Length_Respiratory_System
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Respiratory_System(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Respiratory_System(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Thoracic,1);
        for j = 1:Array_Length_Thoracic
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Thoracic(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Thoracic(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
    Store_Index = Store_Index + 1;
    
    for i = 1:size(Gene_Comparison_Human,1)
        Search_Term = Gene_Comparison_Human(i,2); % Pick a search term from the list
        
        Indices_Human = zeros(Array_Length_Urinary,1);
        for j = 1:Array_Length_Urinary
            Indices_Human(j) = strcmp(Pared_MONARCH_Disease_Urinary(j).Gene_Name, Search_Term); % Find where the term is
        end
        Indices_Human = logical(Indices_Human);             % Convert to logical indices
        
        Disease_Results_Human = {};
        Disease_Results_Human = Pared_MONARCH_Disease_Urinary(Indices_Human);
        Disease_Results_Human = struct2cell(Disease_Results_Human);
        Disease_Results_Human = Disease_Results_Human(4,1,:);
        
        squeeze(Disease_Results_Human);
        size(Disease_Results_Human,3);
        
        Gene_Comparison_Human{i,Store_Index} = size(Disease_Results_Human,3);
    end
 
% -------------------------------------------------------------------------     
    for j = 3:Store_Index
        Summary_Table(x,j-2) = sum([Gene_Comparison_Human{:,j}]);
    end
end

Gene_Comparison_Human
Summary_Table


%% Remove extraneous files, save the final final, then convert to text
clearvars -except Summary_Table
save('Home Directory/Summary_Table.mat')
dlmwrite('Monarch Disease Heat Map Summary.txt',Summary_Table,'delimiter','\t')

