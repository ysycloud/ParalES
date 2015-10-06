
pert_type = 'trt_cp';
cell_id_set={'VCAP','PC3','A375','HA1E','A549','MCF7','HT29','HEPG2','HCC515'};
%pert_desc_set = {'Hydroxytacrine maleate (R,S)','DONEPEZIL HYDROCHLORIDE','donepezil','Rivastigmine','Galanthamine hydrobromide','galanthamine','memantine','Memantine Hydrochloride'};
file_name = 'data_drug265998_level3.txt';  %out file name 

%ds = parse_gctx('/home/yinxy/data/gex_epsilon_n1429794x978.gctx');
ds = parse_gctx('../data/modzs_n272x978.gctx');
[m,n]=size(ds.mat);
getSample 		%get Sample in L1000 attribution before

probe=1:m;            %probe index
probe=probe'; 
   
fid = fopen(file_name, 'w');  %file point source
fprintf(fid,'%10g\t%10g\n', count,m);
tic
%pre-sort
o=ones(m,2);
for i = 1:count
    o = [mat(:,i),probe];
    o = sortrows(o,1);
    mat(:,i)=o(:,2);
    for j = 1:m-1
         fprintf(fid,'%5g\t',o(j,2));
     end
     fprintf(fid,'%5g\n',o(m,2));
end
toc
fclose(fid);

