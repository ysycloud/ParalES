threadnum=50;    %coreNum
sig=50;           		%signature length
per_type1 = 'trt_cp';		%sourcr attribute
per_type2 = 'trt_cp';		%mat attribute
cell_id_set={'VCAP','PC3','A375','HA1E','A549','MCF7','HT29','HEPG2','HCC515'};  %cell_id set

ds = parse_gctx('../data/modzs_n272x978.gctx');
[m,n]=size(ds.mat);
getSampleforMat

probe=1:m;            %probe index
probe=probe';        

%pre-sort
o=ones(m,2);
for i = 1:count1
    o = [source(:,i),probe];
    o = sortrows(o,1);
    source(:,i)=o(:,2);  
end
for i = 1:count2
    o = [mat(:,i),probe];
    o = sortrows(o,1);
    mat(:,i)=o(:,2);  
end

%openµÄmatlabpool
if matlabpool('size')<=0 
    matlabpool('open','local',threadnum);
else
    disp('Already initialized');
end

tic
ES_score=ESScore(source,mat,sig,threadnum);
toc
%save ES.mat ES_score;

matlabpool close;