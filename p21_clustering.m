%% cluster p21 dynamics

% p21_cleanmean and p53_cleanmean are the concentration of p21 and p53 in
% cells
p21 = p21_cleanmean{1}(druginsertion:end,:); % druginsertion represents the time point when we irradiate cells
data_smooth = sgolayfilt(p21,3,11);  % smooth data / remove noise

% cluster p21 dynamics using shape-based distance and binary tree
% see dynClassifier and kshpae functions for details
[subgroup_ID centroid_ID distM] = dynClassifier(data_smooth,2,50);

%% Data normalization - a. normalized to basal level
    % define most freqent values before damage induction as basal level of p21.
clear basal_level_p21 p21_norm

for i=1:length(p21_cleanmean)
    p21= p21_cleanmean{i}(1:druginsertion,:)
    sgf = sgolayfilt(p21,3,11);
    basal_level_p21{i} = mode(sgf);
    p21_norm{i} = p21_cleanmean{i}./repmat(basal_level_p21{i},[size(p21_cleanmean{i},1),1]);
end
clear p21 sgf

%% Data normalization - a. normalized to basal level
    % define most freqent values before damage induction as basal level of p53.
clear basal_level_p53 p53_norm

for i=1:length(p53_cleanmean)
    p53 = p53_cleanmean{i}(1:druginsertion,:)
    sgf = sgolayfilt(p53,3,11);
    basal_level_p53{i} = mode(sgf);
    p53_norm{i} = p53_cleanmean{i}./repmat(basal_level_p53{i},[size(p53_cleanmean{i},1),1]);
end
clear p53 sgf

%% plot single cell heatmap of p53 and p21 by clusters

% set xtick and labels for plots
a = sort(unique([fliplr(Xaxis(druginsertion):-5:0),Xaxis(druginsertion):5:max(Xaxis)]));
xtick = a*60/(deltaT)+1;
xaxis_label = a-Xaxis(druginsertion);


figure(1),
for i =1:2
    
    [B I] = sort(subgroup_ID(:,i));
    subplot(2,2,i*2-1)
    imagesc(zscore(p21_norm{1}(:,I))');colorbar
    ax = gca;
    ax.XTick = xtick
    ax.XTickLabel = xaxis_label;
    title('p21 fold-change')
    caxis([-1 2])
    ylabel('cell #');xlabel('time after irradiation(h)')
    
    subplot(2,2,i*2)
    imagesc(zscore(p53_norm{1}(:,I))');colorbar
    caxis([-1 3])
    ax = gca;
    ax.XTick = xtick;
    ax.XTickLabel = xaxis_label;
    title('p53 fold-change')
    ylabel('cell #');xlabel('time after irradiation(h)')
    
end

print(figure(1),'-dpdf', '-noui',['plot/','21.p53_p21_heatmap ','(',date,').pdf'])

%% plot population level of p53 and p21 by clusters

% plot all cells in 1st column
subplot(4,3,4)
plot2p(p53_cleanmean{1},p21_cleanmean{1},Xaxis,druginsertion) 
% plot2p is a data visualization tool to plot 25th percentile, 75th
% percentile and median level of both p53 and p21 in a single plot
title('all cells')

hold on,

% plot 1st level clusters in 2nd column and 2nd level clusters in 3rd
% column
for ii =1:2
        temp_p21=[]; temp_p53=[]
        temp_p21 = p21_cleanmean{1}(:,subgroup_ID(:,1) ==ii);
        temp_p53 = p53_cleanmean{1}(:,subgroup_ID(:,1) ==ii);
        
        for iii=0:2
            clear temp_p21_subgroup temp_p53_subgroup
            if iii ==0
                temp_p21_subgroup = temp_p21;
                temp_p53_subgroup = temp_p53;
            else
                temp_p21_subgroup = p21_cleanmean{1}(:,subgroup_ID(:,2)==ii*10+iii);
                temp_p53_subgroup = p53_cleanmean{1}(:,subgroup_ID(:,2)==ii*10+iii);
            end
        
        subplot(4,3,ii*6-4+iii*iii)
        hold on,
        title(['subgroup: ',num2str(ii*10+iii)])
        plot2p(temp_p53_subgroup,temp_p21_subgroup,Xaxis,druginsertion)
        end
end

print(gcf,'-dpdf', '-noui',['plot/','1.kshapse_clusters_population_level','(',date,').pdf'])
