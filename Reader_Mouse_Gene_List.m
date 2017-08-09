%% Import data from spreadsheet

[~, ~, Mouse_Gene_List] = xlsread('Home Directory/90 Gene List 070517.xlsx','Sheet1');
Mouse_Gene_List(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),Mouse_Gene_List)) = {''};

save('Home Directory/Mouse_Gene_List.mat')

clearvars Mouse_Gene_List