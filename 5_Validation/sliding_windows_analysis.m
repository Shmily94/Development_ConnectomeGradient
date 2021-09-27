%% sliding window analysis

load('...\sliding_win_ana\order.mat'); % sort subjects according to age;
load('...\table_model_resort.mat');
realigned_gradient_sort_age=zeros(17673,13,491);
age=table_model_resort.age;
for i=1:491
    index=order(i);
    realigned_gradient_sort_age(:,:,i)=realigned_gradient.realigned(:,:,index);
    realigned_age(i,:)=age(index,:);
    realigned_explan(i,:)=gradient_indi_explan_resort(index,:);
end
save('realigned_age','realigned_age')
save('realigned_explan','realigned_explan')
save('realigned_gradient_sort_age','realigned_gradient_sort_age');

%divide windows
win_length=30;
win_step=5;
for i= 1:93
  gradient_window = realigned_gradient_sort_age(:,1:2,1+(i-1)*win_step:win_length+(i-1)*win_step);
  save(['gradient_window',mat2str(i)],'gradient_window');
  gradient1_window(:,i)=mean(gradient_window(:,1,:),3);
  gradient2_window(:,i)=mean(gradient_window(:,2,:),3);
  age_window = realigned_age(1+(i-1)*win_step:win_length+(i-1)*win_step,:);
  mean_age_window(i,:)=mean(age_window);
  explan_window = realigned_explan(1+(i-1)*win_step:win_length+(i-1)*win_step,:);
  mean_explan_window(i,:)=mean(explan_window);
  
end
save('gradient1_window_l30_s5','gradient1_window');
save('gradient2_window_l30_s5','gradient2_window');
save('mean_age_window','mean_age_window');
save('mean_explan_window','mean_explan_window');





