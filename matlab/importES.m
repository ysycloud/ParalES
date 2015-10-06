clear ;
data_fname = '../data/ESResult.dat' ; % file name

file_id = fopen(data_fname, 'rb');
[mn,ele_count] = fread(file_id, 2, 'int32'); % read sample number
row = mn(1);
col = mn(2);

ES_Score = [1:col];

while feof(file_id) == 0

    [row_array, ele_count] = fread(file_id, col, 'float32') ;  %read one row£¬but is col*1
	
    if ele_count < col % elecount < col represent to the file end
        break ;
    else
        % transpose to 1*col
        row_array = row_array'  ;
        % append row_array to raw_data
        ES_Score = [ES_Score; row_array] ;
    end
end
% get off the first line 
ES_Score(1,:)=[] ;
% close file
fclose(file_id);
% delete other usless vars
clear data_fname file_id fid mn ele_count i m n row_array ans;