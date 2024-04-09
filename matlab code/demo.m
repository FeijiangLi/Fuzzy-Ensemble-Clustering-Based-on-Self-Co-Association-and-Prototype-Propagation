load seed.mat
gt=data(:,end);
H=31;
k=length(unique(gt));
data_feature=data(:,1:end-1);
data_feature=predata(data_feature);

[clusters,allk] =creat_clusters_fix_fuzzy(data_feature,H,k);
[cl1,result,initial,pdf] = SCPP_a(clusters,allk,k);
[ac1,ARI1,NMI1]=evaluate2(cl1,gt,k)

[cl2,result,initial] = SCPP_v(clusters,allk,k);
[ac2,ARI2,NMI2]=evaluate2(cl2,gt,k)