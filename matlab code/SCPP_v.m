function [cl,result,initial] = SCPP_v(clusters,allk,k)
n=size(clusters,1);
l=length(allk);

[simf ] = sim_fuzzy(clusters,allk);
sim=sum(simf,3)./l;
% sim=clusters*clusters'./l;
dens=diag(sim);

[lden] = local_dens(1-sim,dens);

use_lden=find(dens>=mean(dens));
[~,use_initial]=sort(lden(use_lden),'descend');
initial=use_lden(use_initial);

in_k=floor(sqrt(n));

if in_k<k
    in_k=k;
end

initial=initial(1:in_k);


[result] = link(simf,initial,in_k);

if in_k==k
    cl=result;
else
    simij=zeros(1,(in_k-1)*(in_k-2)/2+(in_k-1));
    for i=1:in_k
        locat_i=result==i;
        for j=i+1:in_k
            locat_j=result==j;
            sim_ij=sim(locat_i,locat_j);
            ij=sum(sum(sim_ij))/(sum(locat_i)*sum(locat_j)); 
            simij((in_k-1)*(i-1)-((i-1)*(i-2)/2)+j-i)=ij;
        end
    end
    simij=simij./max(simij);
    dis=1-simij;
    Z = linkage(dis,'single');
    clu=cluster(Z,k);
    CC=zeros(n,1);

    for i=1:in_k
        biaoji=clu(i);
        l= result==i;
        CC(l,:)=biaoji;
        
    end
cl=CC;
end

function [lden] = local_dens(dis,dens)
n=length(dens);
lden=zeros(1,n);
maxdis=max(max(dis))+1;

lo_maxdens=find(dens==max(dens));
lden(lo_maxdens(1))=maxdis+1;
if length(lo_maxdens)>1
    dens(lo_maxdens(2:end))=dens(lo_maxdens(2:end))-0.0001;
end

for i=[1:lo_maxdens(1)-1, lo_maxdens(1)+1:n]
    now_dis=dis(i,:);
    now_dis(i)=maxdis+1;
    lo_deni=dens==dens(i);
    
    if sum(lo_deni)>1
        dens(lo_deni)=dens(lo_deni)-0.000001;
        dens(i)=dens(i)+0.000001;
    end
    
    lo= dens>dens(i);
    no_dis=min(now_dis(lo));
    lden(i)=no_dis(1);
end


function [cl] = link(simf,initial,k)
[n,~,l]=size(simf);
cl=zeros(n,1);
sqrtns=ceil(sqrt(n));

for i=1:k
    cl(initial(i))=i;

end
while sum(cl==0)>0
    liu=find(cl==0);
    labeled=n-length(liu);
    sim_in=zeros(length(liu),length(labeled),l);
    for i=1:k
        lo_label_i=cl==i;
        liusimi=simf(liu,lo_label_i,:);
        sim_in(:,i,:)=max(liusimi,[],2);
    end
    
    [m]=max(sim_in,[],2);
    voting=sim_in==m;
    voting=sum(voting,3);
    re_sim_in=sort(voting,2,'descend');
    cha_sim_in=abs(re_sim_in(:,1)-re_sim_in(:,2));
    
    if sum(sum(cha_sim_in))==0
        now_label_lo=1:1:length(cha_sim_in);
        cl(liu(now_label_lo))=k+1;
    else
    
%         th=graythresh(cha_sim_in);
%         now_label_lo=find(cha_sim_in>th);

        [~,lo]=sort(cha_sim_in,'descend');
        now_label_lo=lo(1:min(sqrtns,length(lo)));

        if isempty(now_label_lo)
            now_label_lo=1:1:length(cha_sim_in);
        end
        [~,label]=max(voting(now_label_lo,:),[],2);
        cl(liu(now_label_lo))=label;

    end
end

function [simf ] = sim_fuzzy(cls,allk)
n=size(cls,1);
l=length(allk);
simf=zeros(n,n,l);
sum_k= cumsum(allk);
for i=1:l
    if i==1
        n_cls=cls(:,1:sum_k(i));
    else
        n_cls=cls(:,sum_k(i-1)+1:sum_k(i));
    end
    simf(:,:,i)=n_cls*n_cls';
    
end

