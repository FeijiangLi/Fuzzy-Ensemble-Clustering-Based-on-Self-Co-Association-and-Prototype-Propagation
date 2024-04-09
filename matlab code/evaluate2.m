function [ac,ARI,NMI]=evaluate2(cluster,Ak,k)
%-----------------------------------------------------------------
%     输入：
%          cluster：原始类标签
%          Ak：聚类结果
%          k：类个数
%     输出：
%          ac：精度
%          pre：纯度
%          re：召回率
%          ARI：调整兰德指数
%------------------------------------------------------------------
n=length(cluster);
cp=crosstab(cluster,Ak);
cp(k+1,:)=sum(cp);
cp(:,k+1)=sum(cp,2);
%for j=1:k
%    [a,b]=size(find(A1==j));
%    p(j)=a;
%end
%for j=1:k
%    [a,b]=size(find(A2==j));
%    c(j)=a;
%end
%-----------------------------------------------------------

cpp=cp(1:k,1:k);

a=0;
share_num=zeros(k,k);
for i=1:k  %计算精度
    a=a+max(max(cpp));
    [max_hang,max_lie]=find(cpp==max(max(cpp)));
    share_num(max_hang(1),max_lie(1))=max(max(cpp));
    cpp(max_hang(1),:)=-1;
    cpp(:,max_lie(1))=-1;
end
ac=a/n;

% ac=sum(max(cpp))/n;

%---------------------------------------------------------------
%-----------------------------------------------
%计算NMI

cm=cp(1:end-1,1:end-1);
locat=cm==0;
cp1=cp(1:end-1,end);
cp2=cp(end,1:end-1);
cp12=cp1*cp2;
NMI_cm=cm.*log((n.*cm)./cp12);
NMI_cm(locat)=0;
NMI_up=sum(sum(NMI_cm));
NMI_cp1=sum(cp1.*log(cp1./n));
NMI_cp2=sum(cp2.*log(cp2./n));
NMI_down=sqrt(NMI_cp1*NMI_cp2);
NMI=NMI_up/NMI_down;


% %计算RI
% cp=cp(1:end-1,1:end-1);
% cpO=sum(cp);
% cOp=sum(cp,2);
% 
% cp2=cp.^2;
% cpO2=cpO.^2;
% cOp2=cOp.^2;
% 
% 
% n11=(1/2)*sum(sum(cp.*(cp-1)));
% n00=(1/2)*(n^2+sum(sum(cp2))-(sum(cpO2)+sum(cOp2)));
% n01=(1/2)*(sum(cpO2)-sum(sum(cp2)));
% n10=(1/2)*(sum(cOp2)-sum(sum(cp2)));
% 
% ARI=(n11+n00)/(n11+n00+n10+n01);

%%计算ARI
cp1=cp(1:k,1:k);
r0=sum(sum((cp1.*(cp1-1))./2));
cp2=cp(1:k,1+k)';
r1=sum((cp2.*(cp2-1))./2);
cp3=cp(1+k,1:k);
r2=sum((cp3.*(cp3-1))./2);
r3=(2*r1*r2)/(n*(n-1));
ARI=(r0-r3)/(0.5*(r1+r2)-r3);

end