function [ES_Score] = ESScore(source,mat,sig,corenum)
[m,ns]=size(source);
[~,n]=size(mat);

ES_Score=ones(ns,n);

if ns < corenum    %Fine-grained parallelism
    for i=1:ns
         parfor j = 1:n    
	     ES_Score(i,j)=ESquick(source(:,i) ,mat(:,j),sig);
         end
    end
else         %Coarse-grained parallelism
    parfor i=1:ns
        for j = 1:n    
	    ES_Score(i,j)=ESquick(source(:,i) ,mat(:,j),sig);
        end
    end
end
