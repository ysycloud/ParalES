%get cid list in info
tic
load info.mat;
info = regexp(data,'\t','split');
toc
clear data;

[im,in] = size(info);
[cellm,celln] = size(cell_id_set);
ccid = 0;
untrt_ccid = 0;
tic
for i=2:im
	if isequal(info{i,1}{1,24},pert_type)&&~isequal(info{i,1}{1,3},'-666')&&ismember(info{i,1}{1,4},cell_id_set)   %&&ismember(info{i,1}{1,3},pert_desc_set)
		ccid = ccid + 1;   %find the cids and mark the cell id vector
		info_cid{ccid} = info{i,1}{1,1};
		for j=1:celln
			if ismember(info{i,1}{1,4},cell_id_set(j))
				trtcellflag(ccid) = j; 
				break;
			end		
		end
	end
	if isequal(info{i,1}{1,24},'ctl_untrt')&&~isequal(info{i,1}{1,3},'-666')&&ismember(info{i,1}{1,4},cell_id_set)
        untrt_ccid = untrt_ccid + 1;   %find the cids and mark the cell id vector
		untrt_cid{untrt_ccid} =  info{i,1}{1,1};
		for j=1:celln
			if ismember(info{i,1}{1,4},cell_id_set(j))
				untrtcellflag(untrt_ccid) = j; 
				break;
			end
		end		
	end
end
toc


%compute the averge untrt vector
tic
[untrt_local,untrt_index] = ismember(ds.cid,untrt_cid);
untrtnum = zeros(1,celln);
untrtMat = zeros(m,celln);
for i=1:n
	if untrt_local(i)==1
		current_cellid = untrtcellflag(untrt_index(i));
		untrtnum(1,current_cellid) = untrtnum(1,current_cellid)+1;
		untrtMat(:,current_cellid) = untrtMat(:,current_cellid)+ ds.mat(:,i);
		if sum(untrtnum)>=untrt_ccid
			break;
		end
	end
end
for i=1:celln
	untrtMat(:,i) = untrtMat(:,i)/untrtnum(1,i);  %average
end
toc


%get the sample by cid
tic
[trt_local,trt_index]=ismember(ds.cid,info_cid); 
count = 0;
for i=1:n
	if trt_local(i)==1
		current_cellid = trtcellflag(trt_index(i));
		count = count+1;
		mat(:,count) = ds.mat(:,i)./untrtMat(:,current_cellid);  %°´infoÖĞË³ĞòÅÅ
		cid(count) = ds.cid(i);
		if count>=ccid
			break;
		end
	end
end
toc

clear info 
save cid_drug265998.mat cid
