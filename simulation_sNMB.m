
clear;
clc;
close all;


e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.5;

p=zeros(1,25);

p(1)=-0.5;
p(2)=-0.475;
p(3)=-0.45;
p(4)=-0.4;
p(5)=-0.375;
p(6)=-0.35;
p(7)=-0.325;
p(8)=-0.3;
p(9)=-0.275;
p(10)=-0.25;
p(11)=-0.225;
p(12)=-0.2;
p(13)=-0.175;
p(14)=-0.15;
p(15)=-0.125;
p(16)=-0.1;
p(17)=-0.08;
p(18)=-0.05;
p(19)=-0.02;
p(20)=-0.001;
p(21)=0.02;
p(22)=0.05;
p(23)=0.08;
p(24)=0.1;
p(25)=0.15;

adjacent_network= [1 2 3 4 5 6;...
    2 1 3 4 5 6;...
    3 1 2 0 0 0;...
    4 1 2 0 0 0;...
    5 1 2 6 0 0;...
    6 1 2 5 7 8;...
    7 6 8 0 0 0;...
    8 6 7 0 0 0];


D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
T=1;
TT=100;

reference_sample_num=8;
CC=zeros(8,reference_sample_num);
q(1)=0.96^(1/abs(p(1)));
 E=[-2/5*q(1) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
J=D*E*inv(D);
for k=1:reference_sample_num
    for i=1:N-1;
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8;
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
    CC(:,k)=X(:,2000);
end
pre_TC=CC;
total_node_num=8;

 
for s=1:TT
    for l=2:25
        q(l)=0.96^(1/abs(p(l)));
        E=[-2/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        for i=1:N-1
            ts(i+1)=ts(i)+delta_t;
            eJ=e^(J*delta_t);
            for jj=1:8
                X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
            end
        end
        TC=X(:,2000);
        for na=1:total_node_num
             edges_list=[];
             center=adjacent_network(na,1);
             e=0;
             index_no_zeros=find(adjacent_network(na,:)~=0);
             for n=2:length(index_no_zeros)
                 nei=adjacent_network(na,n);
                 e=e+1;
                edges_list(e,:)=[center nei];
             end
            clear delt_pcc;
            clear delt_sd;
            if e<2
                Land_entropy(na,l-1,s)=0;
                continue;
            end
            for i=1:e
                curr_pcc=abs(corr(pre_TC(edges_list(i,1),:)',...
                    pre_TC(edges_list(i,2),:)'));
                
                temp_add_onecase1=[pre_TC(edges_list(i,1),:),...
                    TC(edges_list(i,1))];
                temp_add_onecase2=[pre_TC(edges_list(i,2),:),...
                    TC(edges_list(i,2))];
                curr_pcc_add_onecase=abs(corr(temp_add_onecase1',temp_add_onecase2'));
                delt_pcc(i)=abs(curr_pcc_add_onecase-curr_pcc);
                delt_sd(i)=abs(std([pre_TC(edges_list(i,2),:),TC(edges_list(i,2))])-std(pre_TC(edges_list(i,2),:)));
            end
            delt_pcc_index=find(delt_pcc~=0);
            if length(delt_pcc_index)==0 || length(delt_pcc_index)==1
                Land_entropy(na,l-1,s)=0;
                continue;
            end
            delt_pcc_mean=mean(delt_pcc);
            delt_sd=[delt_sd,abs(std([pre_TC(center,:),TC(center)])-std(pre_TC(center,:)))];
            delt_sd_mean=mean(delt_sd);
            Land_entropy(na,l-1,s)=delt_pcc_mean*delt_sd_mean;
        end
    end
    s
end
%%
Land_entropy_size=size(Land_entropy);

for s=1:Land_entropy_size(3)
    r=Land_entropy(:,:,s);
    r_size=size(r);
    for t=1:r_size(2)
        a=sort(r(:,t),'descend');
        result(s,t)=mean(a(1:8));
    end
end
sDNBscore=mean(result);


p(2:25)=1:24;
plot(p(2:25),sDNBscore,'r-*','LineWidth',2.5);
xlabel('parameter p');
ylabel('pagerank scoer');
%ylim([-0.01,0.4])
%set(gca,'ytick',[0:0.1:0.4]);
% ylim([0 0.05]);
set(gca,'XTick',[p(2:2:24)]);
B={ '-0.45'  '-0.4','-0.35 ','-0.3','-0.25','-0.2','-0.15','-0.1','-0.05','-0.0001','0.05','0.1'};
set(gca,'XTickLabel',B);










