%% spatial autocorrelation correlation

cd ...\GO_Results
for i = 1:10000
    id = num2str(i-1);
    cd(['...\GO_Results\',id])
    if exist('GOCOMPONENT.xlsx','file') ~=0
    [num,txt,raw]= xlsread('GOCOMPONENT.xlsx');
    [r,c]=size(num);
    surr_GO=cell(r,1);
    for j=1:r
    surr_GO{j}=raw{j+1,2};
    end
    %nan
    temp1=cellfun(@(x) length(find(isnan(x))),surr_GO,'uniformoutput',false);

    temp2=cell2mat(temp1);

    index=find(temp2==1); 
    surr_GO(index)=[]; %#ok<FNDSB>
    [C,ia,ib]=intersect(surr_GO,real_GO);
    for k=1:length(C)
    if isempty(C{k})==0
        surr_real_GO{i}.name=C;
        surr_real_GO{i}.surr_p=num(ia,1);
        surr_real_GO{i}.surr_q=num(ia,2);
        
    end
    end
    end

end


meed=~cellfun('isempty',surr_real_GO);
index2=find(meed==1);
surr_GO_component=surr_real_GO(index2);

n1=0;
n2=0;
n3=0;
n4=0;
n5=0;
for x=1:length(surr_GO_component)
    for y=1:length(surr_GO_component{x}.name)
    str1=surr_GO_component{x}.name(y,1);
    if strcmp(str1,'axon part')==1
        n4=n4+1;
        n4_p(x,1)=surr_GO_component{x}.surr_p(y,1);
    end
    if strcmp(str1,'synapse part')==1
        n5=n5+1;
        n5_p(x,1)=surr_GO_component{x}.surr_p(y,1);

    end
%     if strcmp(str1,'regulation of trans-synaptic signaling')==1
%     n3=n3+1;
%     n3_p(x,1)=surr_GO_component{x}.surr_p(y,1);
% 
%     end
    end 
end
find(n4_p<0.00000491 & n4_p>0);
find(n5_p<0.00000778 & n5_p>0);
    
 
    
