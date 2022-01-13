
clear;
clc;
close all;
[profile,pipi]= xlsread('READ_profile.xlsx');


normal=profile(:,1:11);
case_mprofile=profile(:,12:177);
fid=fopen('READ_tmp_network_adj_edges_all.txt');
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
total_node_num=j;


reference_sample_num=11;
patients_num=[34,53,53,26];
case_mprofile=fillmissing(case_mprofile,'constant',0);
es=0.00000001;

tempcase(:,1,1:patients_num(1))=case_mprofile(:,1:34);     % Stage I 
tempcase(:,2,1:patients_num(2))=case_mprofile(:,35:87);   % Stage II
tempcase(:,3,1:patients_num(3))=case_mprofile(:,88:140);   % Stage III
tempcase(:,4,1:patients_num(4))=case_mprofile(:,141:166);   % Stage IV




psize=size(tempcase);
Land_entropy=zeros(total_node_num,53,4);
stage=4;

for l=1:psize(2)
    for  s=1:patients_num(l)
         for na=1:total_node_num
             edges_list=[];
             center=adjacent_network{na}{1};
             e=0;
             for n=2:length(adjacent_network{na})
                 nei=adjacent_network{na}{n};
                 e=e+1;
                edges_list(e,:)=[str2num(center) str2num(nei)];
             end
                         % psize(3)
            clear delt_sd;
            clear delt_pcc;
            if e<2
                Land_entropy(na,s,l)=0;
                continue;
            end
            for i=1:e
                curr_pcc=abs(corr(normal(edges_list(i,1),:)',...
                    normal(edges_list(i,2),:)'));
                
                temp_add_onecase1=[normal(edges_list(i,1),:),...
                    reshape(tempcase(edges_list(i,1),l,s),1,1)];
                temp_add_onecase2=[normal(edges_list(i,2),:),...
                    reshape(tempcase(edges_list(i,2),l,s),1,1)];
                curr_pcc_add_onecase=abs(corr(temp_add_onecase1',temp_add_onecase2'));
                delt_pcc(i)=abs(curr_pcc_add_onecase-curr_pcc);
                delt_sd(i)=abs(std([normal(edges_list(i,2),:),reshape(tempcase(edges_list(i,2),l,s),1,1)])-std(normal(edges_list(i,2),:)));
            end
            delt_pcc_index=find(delt_pcc~=0);
            if length(delt_pcc_index)==0 || length(delt_pcc_index)==1
                Land_entropy(na,s,l)=0;
                continue;
            end
            delt_pcc_mean=mean(delt_pcc);
            delt_sd=[delt_sd,abs(std([normal(str2num(center),:),reshape(tempcase(str2num(center),l,s),1,1)])-std(normal(str2num(center),:)))];
            delt_sd_mean=mean(delt_sd);
            Land_entropy(na,s,l)=delt_pcc_mean*delt_sd_mean;
         end
    end
end
    
   
save READ_Land_entropy.mat;

Land_entropy_size=size(Land_entropy);
case_result=zeros(53,4);
for t=1:Land_entropy_size(3)
    for case_num=1:patients_num(t)
        [sort_Land_entrop,idx]=sort(Land_entropy(:,case_num,t),'descend');
        case_result(case_num,t)=sum(sort_Land_entrop(1:300));                    
    end
    result(t)=mean(case_result(1:patients_num(t),t));
    
end
figure;
t=[1 2 3 4];
plot(t,result,'r','LineWidth',3);
set(gca,'XTick',1:4);
B={'I' 'II'  'III' 'IV'};
set(gca,'XTickLabel',B);
xlabel('Stages');
ylabel('Entropy');
title('Average Entropy for READ ');












