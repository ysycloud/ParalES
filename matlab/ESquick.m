function [ES] = ESquick(s1,s2,sig)
s=[s1,s2];  
[m,n] = size(s);
score = ones(n,1);
for i=1:n
    if i==1
       up = s1(1:sig,1);
	down = s1((m-sig+1):m,1);
	rank = s2(:,1);
    else
	up = s2(1:sig,1);
	down = s2((m-sig+1):m,1);
	rank = s1(:,1);
    end
     
    isgsUp = ismember(rank,up);
	isgsDown = ismember(rank,down);	

    % compute ESUP
	score_hit	= cumsum(isgsUp);
	score_hit	= score_hit/score_hit(end);
	score_miss  = cumsum(1-isgsUp);
	score_miss  = score_miss/score_miss(end);
    [~,t] = max(abs(score_hit - score_miss));
    ESUP = score_hit(t)-score_miss(t);
       
	% compute ESDOWN
	score_hit	= cumsum(isgsDown);
	score_hit	= score_hit/score_hit(end);
	score_miss  = cumsum(1-isgsDown);
	score_miss  = score_miss/score_miss(end);
    [~,t] = max(abs(score_hit - score_miss));
    ESDOWN = score_hit(t)-score_miss(t);

	%Average ES
    score(i,1)=(ESUP - ESDOWN)/2;
end
ES=(score(1,1)+score(2,1))/2;
return
end