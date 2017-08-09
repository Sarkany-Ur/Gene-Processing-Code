function [Pared_Data_Set_Out] = Data_Paring(Data_Set_In)

Data_Set_In = table2cell(Data_Set_In);        % Convert from table to cell
Array_Length = size(Data_Set_In,1);           % Find the length of the data set

Pared_Data_Set_Out(Array_Length) = struct();  % Use a structure to store

for i = 1:Array_Length
    Pared_Data_Set_Out(i).NCBI_Gene_ID = str2double(char(Data_Set_In{i,1}(1,10:end)));
    Pared_Data_Set_Out(i).Gene_Name    = char(Data_Set_In{i,2});
    Pared_Data_Set_Out(i).HPO_Number   = str2double(char(Data_Set_In{i,3}(1,4:end)));
    Pared_Data_Set_Out(i).HPO_Terms    = char(Data_Set_In{i,4});
end


end
