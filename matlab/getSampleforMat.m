
%get cid list in info
tic
load info.mat;
info = regexp(data,'\t','split');
clear data£»
toc
[im,in] = size(info);
[cellm,celln] = size(cell_id_set);
ccid1 = 0;
ccid2 = 0;
untrt_ccid = 0;
tic
for i=2:im
	if isequal(info{i,1}{1,24},pert_type1)&&~isequal(info{i,1}{1,3},'-666')&&ismember(info{i,1}{1,4},cell_id_set)
		ccid1 = ccid1+1;
		info_cid1{ccid1,1} = info{i,1}{1,1};
		for j=1:celln
			if ismember(info{i,1}{1,4},cell_id_set(j))
				trtcellflag1(ccid1) = j; 
				break;
			end		
		end		
	end
	if isequal(info{i,1}{1,24},pert_type2)&&~isequal(info{i,1}{1,3},'-666')&&ismember(info{i,1}{1,4},cell_id_set)
		ccid2 = ccid2+1;
		info_cid2{ccid2,1} = info{i,1}{1,1};	
		for j=1:celln
			if ismember(info{i,1}{1,4},cell_id_set(j))
				trtcellflag2(ccid2) = j; 
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


%get the sample  by cid
tic
[trt_local1,trt_index1]=ismember(ds.cid,info_cid1);
[trt_local2,trt_index2]=ismember(ds.cid,info_cid2); 
count1 = 0;
count2 = 0;
for i=1:n
	if trt_local1(i)==1
		current_cellid = trtcellflag(trt_index1(i));
		count1 = count1+1;
		source(:,count1) = ds.mat(:,i)./untrtMat(:,current_cellid);  %°´infoÖĞË³ĞòÅÅ
		cid1(count1) = ds.cid(i);
	end
	if trt_local2(i)==1
		current_cellid = trtcellflag(trt_index2(i));
		count2 = count2+1;
		mat(:,count2) = ds.mat(:,i)./untrtMat(:,current_cellid);  %°´infoÖĞË³ĞòÅÅ
		cid2(count2) = ds.cid(i);
	end
	if count1>=ccid1&&count2>=ccid2
		break;
	end
end
toc

clear info
save cid_mat.mat cid1 cid2
