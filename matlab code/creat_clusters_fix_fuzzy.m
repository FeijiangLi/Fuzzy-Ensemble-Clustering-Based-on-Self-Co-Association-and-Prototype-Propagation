function [clusters,allk] =creat_clusters_fix_fuzzy(data,H,k)
[n,~]=size(data);

maxk=floor(sqrt(n));
mk=min(maxk,50);
k=max(k,mk);

% k=2*k;
% k=3;
% clusters=zeros(n,k*H);
% clusters=zeros(n,k,H);
% ss=rand(1,H);
% 
% allU=(rand(1,H)*(10-1))+1;
% alld=round(rand(1,H)*(10-2))+2;
allk=k.*ones(1,H);

clusters=zeros(n,k,H);
for i=1:H
%     rng(ss(i)*1e8) ;
%     U=allU(i);
%     [~,label]=fcm(data,allk(i),[alld(i),200,NaN,0]);
    [~,label]=fcm(data,k,[2,200,NaN,0]);
%     [~,label]=fcm(data,allk(i));
    
%     clusters(:,k*(i-1)+1:k*i)=label';
    clusters(:,:,i)=label';
%     [~,label]=fcm(data,k);
%     clusters(:,:,i)=label';   
    
end
clusters=reshape(clusters,n,k*H);
end