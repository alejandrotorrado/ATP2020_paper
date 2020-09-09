% MDdropFigureCode
%
% Juliet Bottorff - 2019
%
%
% 
% Code for Plots for Alejandro: Figure 6
%
%
%
%
%%%%% Total percentages for each cell in each state during decline %%%%%
%%%%% AND CORRELATIONS of percent in light and percent in wake and total change in FR

 colors = ['k', 'm', 'b', 'g', 'r', 'c', 'y', 'w'];

 all_MD=[down_states.MD];
 
 anim_names = unique({animal_structure.animal});
 anim_ids = [animal_structure.identifier];
 all_percents=[];
 control_percents=[];

 %names={'lights on', 'lights off', 'REM', 'NREM', 'AW', 'QW'};
 names = {'REM', 'NREM', 'AW', 'QW'};
 c = 1:length(names);

 total_changes=[];
 tot_change_byanim = struct; % deprived only
 persleep_byanim = struct; % deprived only
 total_changes_control=[];
 
 for i=1:length(all_MD)
     id = down_states(i).identifier;
     id = find(anim_ids==id);
     
         percents=down_states(i).totalpercents;

         if down_states(i).MD==1
            all_percents(end+1,:)=percents;
            total_changes(end+1)=down_states(i).total_diffFR;
       
            this_start = animal_structure(id).decline_bins(1)*(60*15);
            this_start = this_start/3600; % in hours
            this_start = this_start - (2.5*24); % in hours from MD
            if this_start<24
                t = 1;
            else
                t = 2;
            end
            
            this_field = [char(animal_structure(id).animal), '_', num2str(t)];
            
             if ~isfield(tot_change_byanim, this_field)
                    tot_change_byanim.(this_field) = [];
             end
             tot_change_byanim.(this_field)(end+1) = down_states(i).total_diffFR;
             if ~isfield(persleep_byanim, this_field)
                    persleep_byanim.(this_field) = [];
             end
             persleep_byanim.(this_field)(end+1) = percents(1) + percents(2);
          
         else
             control_percents(end+1,:)=percents;
             total_changes_control(end+1)=down_states(i).total_diffFR;
      
         end
  
 end
 
 meanpercents=mean(all_percents);
 stdpercents=std(all_percents);
 sempercents=stdpercents./(sqrt(length(all_percents(:,1))-1));
 meanpercents_control=mean(control_percents);
 stdpercents_control=std(control_percents);
 sempercents_control=stdpercents_control./(sqrt(length(control_percents(:,1))-1));
 
all_percents_mean=[];
all_percents_sem=[];
for m=1:length(meanpercents)
    all_percents_mean(m, :)=[meanpercents(m), meanpercents_control(m)];
    all_percents_sem(m, :)=[sempercents(m), sempercents_control(m)];
end
figure()
b = bar (c, all_percents_mean);
set(b(1), 'FaceColor', 'm', 'EdgeColor', 'none')
set(b(2), 'FaceColor', 'k', 'EdgeColor', 'none');
 ylabel('Fraction time during decline')
 set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 20)

 hold on
 
for n=1:length(c)
    
        plot ([c(n)-0.15, c(n)-0.15], [meanpercents(n)+sempercents(n), meanpercents(n)-sempercents(n)], 'k', 'LineWidth', 2)
        hold on
    
        plot ([c(n)+0.15, c(n)+0.15], [meanpercents_control(n)+sempercents_control(n), meanpercents_control(n)-sempercents_control(n)], 'k', 'LineWidth', 2)
        hold on
   
end


persleep = all_percents(:,1) + all_percents(:,2);
perwake = all_percents(:,3) +  all_percents(:,4);
figure()
% MD
scatter(persleep, total_changes, 'm', 'LineWidth', 3)
set (gca, 'box', 'off', 'fontsize', 15)
hold on
xlabel('Percent decline time in Sleep')
ylabel('Change in FR during decline')
title ('Deprived cells')

[p,s,mu]=polyfit(persleep', total_changes, 1);
x=0:.1:1;
y=x.*p(1) + p(2);
plot(x,y,'k', 'LineWidth', 2)
[rho, pval]=corr(persleep, total_changes');
xlim([0.1, 1.1]);
ylim([-1.1, 1.1]);

% Plot correlation by animal instead of by cell: this is for
% statetimes_MD_v6.m where analyzing one of two time periods for every cell
figure()

anim_sleeps = [];
anim_changes = [];
MarkerShapes = {'o', 'd'};

for j = 1:length(anim_names)
    for k = 1:2
        my_field = [char(anim_names(j)), '_', num2str(k)];
        if isfield(tot_change_byanim, my_field)
            mean_animchange = mean(tot_change_byanim.(my_field));
            anim_changes = [anim_changes, mean_animchange];
            sem_animchange = std(tot_change_byanim.(my_field))/sqrt(length(tot_change_byanim.(my_field))-1);
            mean_animsleep = mean(persleep_byanim.(my_field));
            anim_sleeps = [anim_sleeps, mean_animsleep];
            plot(mean_animsleep, mean_animchange, char(MarkerShapes(k)), 'MarkerEdgeColor', 'k', 'MarkerSize', 12, 'MarkerFaceColor', colors(j))
            hold on
            plot ([mean_animsleep, mean_animsleep], [mean_animchange-sem_animchange, mean_animchange+sem_animchange], 'k')
        end
    end
end

set(gca,'unit','Normalized','position',[.1 .1 .6 .5]);
set (gca, 'box', 'off', 'fontsize', 15)
xlabel('Percent decline time in Sleep')
ylabel('Change in FR during decline')
title ('Deprived cells, by Animal and Drop Start Time ')

[p,s,mu]=polyfit(anim_sleeps, anim_changes, 1);
x=0.05:.1:0.9;
y=x.*p(1) + p(2);
plot(x,y,'k', 'LineWidth', 2)
plot([0,1], [0,0], 'k--', 'LineWidth', 2)
[rho, pval]=corr(anim_sleeps', anim_changes');
xlim([0, 1.1]);
ylim([-1.1, 1.1]);

%%
%%%% State Dense analysis Averaged BY CELL %%%%
all_MD=[down_states.MD];

sleepdense=[];
sleepdense_control=[];
wakedense=[];
wakedense_control=[];
lightdense=[];
lightdense_control=[];
darkdense=[];
darkdense_control=[];
REMdense=[];
REMdense_control=[];
NREMdense=[];
NREMdense_control=[];
AWdense=[];
AWdense_control=[];
QWdense=[];
QWdense_control=[];

%names={'Sleep', 'Wake', 'Light', 'Dark'};
names = {'Sleep', 'Wake'};
c = 1:length(names);

for i=1:length(all_MD)
    if down_states(i).MD==1
        sleepdense=[sleepdense, nanmean(down_states(i).sleepdense_sw)];
        wakedense=[wakedense, nanmean(down_states(i).wakedense_sw)];

    else
        sleepdense_control=[sleepdense_control, nanmean(down_states(i).sleepdense_sw)];
        wakedense_control=[wakedense_control, nanmean(down_states(i).wakedense_sw)];

    end
    
    
end
figure()
meanstates=[nanmean(sleepdense), nanmean(wakedense)];
meanstates_control=[nanmean(sleepdense_control), nanmean(wakedense_control)];

semstates=[nanstd(sleepdense)/sqrt(length(sleepdense)-1), nanstd(wakedense)/sqrt(length(wakedense)-1)];
semstates_control=[nanstd(sleepdense_control)/sqrt(length(sleepdense_control)-1), nanstd(wakedense_control)/sqrt(length(wakedense_control)-1)];

all_states=[];
for j=1:length(c)
    all_states(j,:)=[meanstates(j), meanstates_control(j)];
end
barhandle=bar (c, all_states);
ylabel('Change in FR')
hold on
c_deprived = c-0.15;
c_control = c+0.15;
errorbar(c_deprived, meanstates, semstates,'linestyle','none','capsize',10,'color','k');
errorbar(c_control, meanstates_control, semstates_control,'linestyle','none','capsize',10,'color','k');

%{
for j=1:length(meanstates)
    plot([c(j)-0.15 c(j)-0.15], [meanstates(j)-semstates(j), meanstates(j)+semstates(j)], 'k')
    hold on
    plot([c(j)+0.15 c(j)+0.15], [meanstates_control(j)-semstates_control(j), meanstates_control(j)+semstates_control(j)],'k')

end
%}
hold on
h=cell(1,14);
h{1}='Deprived';
h{2}='Control';
legend (h(1:2))
set (barhandle(1), 'FaceColor', 'm')
set (barhandle(2), 'FaceColor', 'k')
set (gca, 'fontsize', 15, 'xticklabels', names,'box', 'off')
set(gca,'unit','Normalized','position',[.1 .1 .6 .7]);
ylim([-0.2, 0.1])
on=wakedense;
off=sleepdense;
[h,p]=ttest2(on, off)


%%
%%%%%% Plotting start of decline for each cell by animal %%%%%
% need animal_structure (saved in fromATP_DATA folder on server)

% 1 = REM
% 2 = NREM
% 3 = Gen Wake ** NOT USED **
% 4 = Active Wake
% 5 = Quiet Wake

names = unique({animal_structure.animal});
colors = ['k', 'm', 'b', 'g', 'r', 'c', 'y', 'w'];

layer_assign(:,1) = setdiff(1:32,[8 16 24 32]);
% layer_assign(:,2) = [1 4 4 3 3 6 6 1 5 6 3 4 4 3 1 5 5 4 3 1 3 6 5 3 4 4 4 1];
layer_assign(:,2) = [6 2 5 3 4 4 3 2 5 4 4 3 3 6 6 2 5 3 4 4 3 2 5 4 4 3 3 6];

decline_starts = struct;
decline_lengths = struct;
initial_FRs = struct;
layers = struct;
channels = struct;


chunk_1 = []; % drop in first 12 hours after MD (DARK)
chunk_2 = []; % drop in next 12 hours (LIGHT)
chunk_3 = []; % drop in next 12 hours (DARK)
chunk_4 = []; % drop in next 12 hours (LIGHT)
chunk_5 = []; % drop in next 12 hours (DARK)

for i = 1:length(animal_structure)
    if animal_structure(i).MD == 1
        this_start = animal_structure(i).decline_bins(1)*(60*15);
        this_start = this_start/3600; % in hours
        this_start = this_start - (2.5*24); % in hours from MD
        if this_start < 12
            chunk_1(end+1, :) = [animal_structure(i).identifier, animal_structure(i).normFR(2,:)];
        elseif this_start < 24
            chunk_2(end+1, :) = [animal_structure(i).identifier, animal_structure(i).normFR(2,:)];
        elseif this_start < 36
            chunk_3(end+1, :) = [animal_structure(i).identifier, animal_structure(i).normFR(2,:)];
        else 
            chunk_4(end+1, :) = [animal_structure(i).identifier, animal_structure(i).normFR(2,:)]; 
        end
        
        this_end = animal_structure(i).decline_bins(2)*(60*15);
        this_end = this_end/3600;
        this_end = this_end - (2.5*24);
        
        if ~isfield(decline_starts, char(animal_structure(i).animal))
            decline_starts.(char(animal_structure(i).animal)) = [];
        end
        decline_starts.(char(animal_structure(i).animal))(end+1) = this_start;
        if ~isfield(decline_lengths, char(animal_structure(i).animal))
            decline_lengths.(char(animal_structure(i).animal)) = [];
        end
        decline_lengths.(char(animal_structure(i).animal))(end+1) = (this_end - this_start);
        if ~isfield(initial_FRs, char(animal_structure(i).animal))
            initial_FRs.(char(animal_structure(i).animal)) = [];
        end
        initial_FRs.(char(animal_structure(i).animal))(end+1) = animal_structure(i).baseline_FR;
        
        if ~isfield(layers, char(animal_structure(i).animal))
            layers.(char(animal_structure(i).animal)) = [];
        end
        if ~isfield(channels, char(animal_structure(i).animal))
            channels.(char(animal_structure(i).animal)) = [];
        end
        chan = animal_structure(i).channel;
        row = find(layer_assign(:,1)==chan);
        layers.(char(animal_structure(i).animal))(end+1) = layer_assign(row,2);
        channels.(char(animal_structure(i).animal))(end+1) = chan;
        
    end
end

% Plot by animal
figure()
rectangle('Position', [0 0 12 80], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 80], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 80], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 80], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 80], 'FaceColor', [0 0 0 0.15])
count = 1;
for j = 1:length(names)
    if isfield(decline_starts, char(names(j)))
    for jj = 1:length(decline_starts.(char(names(j))))
        plot ( decline_starts.(char(names(j)))(jj), count, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', colors(j))
        count = count+1;
        hold on
    end
    end
end

%plot ([12,12], [0, 80], 'k', 'LineWidth', 3)
%plot ([24, 24], [0, 80], 'k', 'LineWidth', 3)
%plot ([36, 36], [0, 80], 'k', 'LineWidth', 3)
%plot ([48, 48], [0,80], 'k', 'LineWidth', 3)
xlabel ('Time After MD (Hours)')
ylabel ('Cell Number')
title ('Start of FR Drop by Cell')
set (gca, 'FontSize', 15, 'box', 'off')

%{
% Plot organized by initial FR
figure()
rectangle('Position', [0 0 12 20], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 20], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 20], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 20], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 20], 'FaceColor', [0 0 0 0.15])

for j = 1:length(names)
    if isfield(decline_starts, char(names(j)))
    for jj = 1:length(decline_starts.(char(names(j))))
        plot ( decline_starts.(char(names(j)))(jj), initial_FRs.(char(names(j)))(jj), 'ko', 'MarkerSize', 15, 'MarkerFaceColor', colors(j))
      
        hold on
    end
    end
end
xlabel ('Time After MD (Hours)')
ylabel ('Baseline FR (spikes/second)')
title ('Start of FR Drop by Cell, ordered by BASELINE FR')
set (gca, 'FontSize', 15, 'box', 'off')


% Plot organized by LENGTH of decline
figure()
rectangle('Position', [0 0 12 40], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 40], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 40], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 40], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 40], 'FaceColor', [0 0 0 0.15])

for j = 1:length(names)
    if isfield(decline_starts, char(names(j)))
    for jj = 1:length(decline_starts.(char(names(j))))
        plot ( decline_starts.(char(names(j)))(jj), decline_lengths.(char(names(j)))(jj), 'ko', 'MarkerSize', 15, 'MarkerFaceColor', colors(j))
      
        hold on
    end
    end
end
xlabel ('Time After MD (Hours)')
ylabel ('Length of FR Drop (Hours)')
title ('Start of FR Drop by Cell, ordered by LENGTH of drop')
set (gca, 'FontSize', 15, 'box', 'off')


% Plot organized by PUTATIVE layer
figure()
rectangle('Position', [0 0 12 8], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 8], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 8], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 8], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 8], 'FaceColor', [0 0 0 0.15])

for j = 1:length(names)
    if isfield(decline_starts, char(names(j)))
    for jj = 1:length(decline_starts.(char(names(j))))
        plot ( decline_starts.(char(names(j)))(jj), layers.(char(names(j)))(jj), 'ko', 'MarkerSize', 15, 'MarkerFaceColor', colors(j))
      
        hold on
    end
    end
end
xlabel ('Time After MD (Hours)')
ylabel ('Visual Cortex Layer')
title ('Start of FR Drop by Cell, ordered by V1 LAYER')
set (gca, 'FontSize', 15, 'box', 'off')
%}

% Plot organized by CHANNEL (proximity of cells)
figure()
rectangle('Position', [0 0 12 35], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 35], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 35], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 35], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 35], 'FaceColor', [0 0 0 0.15])

for j = 1:length(names)
    if isfield(decline_starts, char(names(j)))
    for jj = 1:length(decline_starts.(char(names(j))))
        plot ( decline_starts.(char(names(j)))(jj), channels.(char(names(j)))(jj), 'ko', 'MarkerSize', 15, 'MarkerFaceColor', colors(j))
      
        hold on
    end
    end
end
xlabel ('Time After MD (Hours)')
ylabel ('Electrode Channel #')
title ('Start of FR Drop by Cell, ordered by electrode')
set (gca, 'FontSize', 15, 'box', 'off')

%{
% Plot MIDDLE OF DECLINE period
figure()
rectangle('Position', [0 0 12 80], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 80], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 80], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 80], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 80], 'FaceColor', [0 0 0 0.15])
count = 1;
for j = 1:length(names)
    if isfield(decline_starts, char(names(j)))
    for jj = 1:length(decline_starts.(char(names(j))))
        midpoint = decline_lengths.(char(names(j)))(jj)/2;
        plot ( (decline_starts.(char(names(j)))(jj) + midpoint), count, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', colors(j))
        count = count + 1;
        hold on
    end
    end
end
xlabel ('Time After MD (Hours)')
ylabel ('Cell Number')
title ('MIDDLE of FR Drop by Cell')
set (gca, 'FontSize', 15, 'box', 'off')
%}

% Plot time course of cells that drop in each 12 hour chunk after MD

mean_chunk1 = nanmean(chunk_1(:, 2:end),1);
sem_chunk1 = nanstd(chunk_1(:, 2:end),1)./sqrt((length(chunk_1(1:end,1))-1));
mean_chunk2 = nanmean(chunk_2(:, 2:end),1);
sem_chunk2 = nanstd(chunk_2(:, 2:end),1)./sqrt((length(chunk_2(1:end,1))-1));
mean_chunk3 = nanmean(chunk_3(:, 2:end),1);
sem_chunk3 = nanstd(chunk_3(:, 2:end),1)./sqrt((length(chunk_3(1:end,1))-1));
mean_chunk4 = nanmean(chunk_4(:, 2:end),1); 
sem_chunk4 = nanstd(chunk_4(:, 2:end),1)./sqrt((length(chunk_4(1:end,1))-1));

first24hours = [chunk_1;chunk_2];
mean_first24hours = nanmean(first24hours(:,2:end),1);
sem_first24hours = nanstd(first24hours(:,2:end),1)./sqrt((length(first24hours(1:end,1))-1));
second24hours = [chunk_3;chunk_4];
mean_second24hours = nanmean(second24hours(:,2:end),1);
sem_second24hours = nanstd(second24hours(:,2:end),1)./sqrt((length(second24hours(1:end,1))-1));

x = animal_structure(1).normFR(1,:).*24; % in hours
%x = x-(2.5*24); % align to MD
%{
figure()
rectangle('Position', [0 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 3], 'FaceColor', [0 0 0 0.15])
shadedErrorBar(x, mean_chunk1, sem_chunk1, 'm')
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that START drop in FIRST 12 HOURS AFTER MD')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')

figure()
rectangle('Position', [0 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 3], 'FaceColor', [0 0 0 0.15])
shadedErrorBar(x, mean_chunk2, sem_chunk2, 'm')
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that START drop in HOURS 12-24 AFTER MD')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [(0), (2.5*24)], 'box', 'off')

figure()
rectangle('Position', [0 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 3], 'FaceColor', [0 0 0 0.15])
shadedErrorBar(x, mean_chunk3, sem_chunk3, 'm')
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that START drop in HOURS 24-36 AFTER MD')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')

figure()
rectangle('Position', [0 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 3], 'FaceColor', [0 0 0 0.15])
shadedErrorBar(x, mean_chunk4, sem_chunk4, 'm')
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that START drop in HOURS 36-48 AFTER MD')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')
%{
figure()
rectangle('Position', [0 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 3], 'FaceColor', [0 0 0 0.15])
shadedErrorBar(x, mean_chunk5, sem_chunk5, 'm')
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that START drop in HOURS 48-60 AFTER MD')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')
%}
%}
figure()

rectangle('Position', [36 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [48 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [60 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [72 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [84 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [96 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [108 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [120 0 12 3], 'FaceColor', [1 1 0 0.15])

shadedErrorBar(x, mean_first24hours, sem_first24hours, 'm')
hold on
plot([0, max(x)], [1,1], 'k--', 'LineWidth', 2)

xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that start drop in FIRST 24 HOURS after MD')
set (gca, 'fontsize', 15)
%set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')
set (gca, 'xlim', [36, 120], 'box', 'off')

set(gca,'unit','Normalized','position',[.1 .1 .6 .5]);
ylim([0,2])

figure()

rectangle('Position', [36 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [48 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [60 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [72 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [84 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [96 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [108 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [120 0 12 3], 'FaceColor', [1 1 0 0.15])

shadedErrorBar(x, mean_second24hours, sem_second24hours, 'm')
hold on
plot([0, max(x)], [1,1], 'k--', 'LineWidth', 2)
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that start drop in SECOND 24 HOURS after MD')
set (gca, 'fontsize', 15)
%set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')
set (gca, 'xlim', [36, 120], 'box', 'off')

set(gca,'unit','Normalized','position',[.1 .1 .6 .5]);
ylim([0,2])

%%
% Checking individual cell dynamics
figure()
rectangle('Position', [0 0 12 3], 'FaceColor', [0 0 0 0.15])
hold on
rectangle('Position', [12 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [24 0 12 3], 'FaceColor', [0 0 0 0.15])
rectangle('Position', [36 0 12 3], 'FaceColor', [1 1 0 0.15])
rectangle('Position', [48 0 12 3], 'FaceColor', [0 0 0 0.15])
plot(x, chunk_4(1,2:end), 'm')
xlabel ('Time (Hours)')
ylabel ('Normalized FR')
title ('Cells that START drop in HOURS 12-24 AFTER MD')
set (gca, 'fontsize', 15)
set (gca, 'xlim', [0, (2.5*24)], 'box', 'off')


%%
%%%%% Plotting example cell declines %%%%
% Need MDCELLS and animal_states for this (get animal_structure from
% running MD_FR_changes_v4.m

numDays = 8.5;
bintime=(12*60*60); %bin time in seconds, for baseline
bintime2=(30*60); %bin time in seconds, for calculating FR
bintime3=(15*60);

%ex_cells = [57, 163, 54, 157, 154, 42, 33, 36, 165, 286, 282, 148, 169, 165];
ex_cells = [165, 282]; % early drop, KH40 and KH67
%ex_cells = [54, 157]; % late drop, KH45 and KH36
for m = 1:length(ex_cells)
    
    n=ex_cells(m); % index of current cell in MDCELLS structure

    
    sorted_time=sort(MDCELLS(n).time); % sort all spiketimes to make sure they're in order
         
        % filter out all spikes from 'off' times
        onTimes=sort(MDCELLS(n).onTime);
        offTimes=sort(MDCELLS(n).offTime);
        
        sorted_time(sorted_time<onTimes(1))=NaN;
        sorted_time(sorted_time>offTimes(end))=NaN;
        
        
        if length(onTimes)>1 
            
            for i=2:length(onTimes)
                current_working_time=sorted_time(sorted_time>offTimes(i-1));
                start_next_off=find(sorted_time==min(current_working_time));
                current_working_time=sorted_time(sorted_time<onTimes(i));
                end_next_off=find(sorted_time==max(current_working_time));
                sorted_time(start_next_off:end_next_off)=NaN;
         
            end
        end

        % spiketimes to use for analysis
        ontime_spikes=sorted_time(~isnan(sorted_time));
experiment_time_start=0;
experiment_time_end= (3600*24*numDays);
edges=experiment_time_start:bintime:experiment_time_end;
[bincounts, binedges]=histcounts(ontime_spikes, edges);
 baseline_FR=sum(bincounts(4:5))/(bintime*2); % 'baseline FR' is averaged FR over 36 hours right before MD
                
                if MDCELLS(n).onTime(1)>36*3600
                    baseline_sp = ontime_spikes(ontime_spikes>=MDCELLS(n).onTime(1));
                    baseline_sp = baseline_sp(baseline_sp<=60*3600);
                    if ((60*3600)-MDCELLS(n).onTime(1))<(3600*5) % make sure baseline time is long enough to mean something!
                    disp (n)
                    disp ('This cell baseline too short!')

                    continue
                    else
                        baseline_FR = length(baseline_sp)/((60*3600)-MDCELLS(n).onTime(1));
                    end
                end
            % Normalize all FRs to baseline 

            edges2=experiment_time_start:bintime2:experiment_time_end; % for finding normalized FR
            edges3=experiment_time_start:bintime3:experiment_time_end; % for sliding window to find decline bins
            [bincounts2, binedges2]=histcounts(ontime_spikes, edges2);
            [bincounts3, binedges3]=histcounts(ontime_spikes, edges3);

            FR_normalized=(bincounts3/bintime3)/baseline_FR;
            for ff=1:length(MDCELLS(n).onTime)
                if MDCELLS(n).offTime(ff)<experiment_time_end
                    if length(MDCELLS(n).onTime)>ff
                        FR_normalized(round(MDCELLS(n).offTime(ff)/bintime3):round(MDCELLS(n).onTime(ff+1)/bintime3))=NaN;
                    else
                        FR_normalized(round(MDCELLS(n).offTime(ff)/bintime3):end)=NaN;
                    end
                end
            end
            FR_normalized(FR_normalized>12)=NaN; % exclude normalized FR above 12
            id = find([animal_structure.identifier]==n);
            db = animal_structure(id).decline_bins;
            %db_secs = db.*bintime3;
            
            %decline_FR = FR_normalized(db(1):db(2));

            xaxis = 0:bintime2:(3600*24*numDays);
            xaxishours = xaxis./(3600);
            %{
            figure()
     
            plot(xaxisdays(1:length(FR_normalized)), FR_normalized, 'k', 'LineWidth', 2)
            xlabel ('Time (Days)')
            ylabel ('Normalized FR')
            title (MDCELLS(n).animal)
            set (gca, 'fontsize', 15)
            set (gca, 'xlim', [0, 8], 'box', 'off')
            hold on
            plot([xaxisdays(db(1)), xaxisdays(db(1))], [0,2.5], 'b', 'LineWidth', 2)
            plot([xaxisdays(db(2)), xaxisdays(db(2))], [0,2.5], 'b', 'LineWidth', 2)
            plot(xaxisdays(db(1):db(2)), decline_FR, 'm', 'LineWidth', 2)
            %}
            
           % myfigname = ['Deprived' num2str(n)];
            clear FR_normalized
            FR_normalized = animal_structure(id).normFR(2,:);
            bl_to_show = mean(FR_normalized(round((1.5*24*3600)/(30*60)):round((2.5*24*3600)/(30*60))))*baseline_FR;
                   % figure ('Name', myfigname)
                    %subplot (2,1,1)
                    %{
                     LD_transition = 0:12:168;
                    LD_transition = LD_transition./24;
                    for j = 1:length(LD_transition)
                        if rem(j,2)
                            rectangle('Position', [LD_transition(j) 0 .5 3], 'FaceColor', [1 1 0 0.1])
                            hold on
                        else
                            rectangle('Position', [LD_transition(j) 0 .5 3], 'FaceColor', [0 0 0 0.1])
                            hold on
                        end
                    end
                    %}
                    %yyaxis left
                    %line([min(xaxisdays), max(xaxisdays)], [1,1], 'Color', 'k')
                    plot(xaxishours(1:length(FR_normalized)), FR_normalized.*baseline_FR, 'k')
                    
                    hold on 
                    %plot ([(db(1)*(15*60))/(3600*24), (db(1)*(15*60))/(3600*24)], [0,max(FR_normalized.*baseline_FR)], 'r', 'LineWidth', 3)
                    plot ([36, 120], [bl_to_show, bl_to_show], 'r--')
                    %clear bl_to_show
                    %{
                    for v=1: length(MDCELLS(n).offTime)
                        plot ([MDCELLS(n).offTime(v)/(3600*24), MDCELLS(n).offTime(v)/(3600*24)], [0,3], 'b')
                    end
                    for w=1: length(MDCELLS(n).onTime)
                        plot ([MDCELLS(n).onTime(w)/(3600*24), MDCELLS(n).onTime(w)/(3600*24)], [0,3], 'y')
                    end
                        %}
                    hold on
                    
                    xlabel('Time from MD (days)')
                    
                    ylabel('Firing Rate (Spikes/Second)')
                    set(gca,'yscale','log');
                    set(gca,'unit','Normalized','position',[.1 .1 .6 .5]);
                    ylim([10e-4, 10e1])

                    title(MDCELLS(n).animal)
                    set (gca, 'fontsize', 15, 'box', 'off')
                    set (gca, 'xlim', [1.5*24,5*24])
                    xticks([36, 48, 60, 72, 84, 96,108,120])
                    yticks([10e-3, 10e-2, 10e-1, 10e0, 10e1, 10e2])
                    %xticks([1.5, 2.5, 3.5, 4.5])
                    %xticklabels({'-1', '0', '1', '2'})
                    
                    %{
                    subplot(2,1,2)
                    MD_start_bintime3=round((3600*24*2.73)/bintime3);
                    MD_end_bintime3=MD_start_bintime3+((3600*24*2)/bintime3);
                    x_md = animal_structure(id).FRslope(1,MD_start_bintime3:MD_end_bintime3);
                    low_slope_indices = find(animal_structure(id).FRslope(2,MD_start_bintime3:MD_end_bintime3)<animal_structure(id).slopeThreshold);
                    slopes_below_threshold = animal_structure(id).FRslope(2,MD_start_bintime3:MD_end_bintime3);
                    slopes_below_threshold = slopes_below_threshold(low_slope_indices);
                    high_slope_indices = find(animal_structure(id).FRslope(2,MD_start_bintime3:MD_end_bintime3)>animal_structure(id).posThreshold);
                    slopes_above_threshold = animal_structure(id).FRslope(2,MD_start_bintime3:MD_end_bintime3);
                    slopes_above_threshold = slopes_above_threshold(high_slope_indices);
                    
                    yyaxis left
                    hold on
                    plot(animal_structure(id).FRslope(1,:), animal_structure(id).FRslope(2,:), 'b')
                
                    plot(x_md(low_slope_indices), slopes_below_threshold, 'xr')
                 
                    plot (x_md(high_slope_indices), slopes_above_threshold, 'xc')
                   
                    plot ([(db(1)*(15*60))/(3600*24), (db(1)*(15*60))/(3600*24)], [0,2.5], 'r', 'LineWidth', 3)
                    
                    xlabel('Time(days)')
                    ylabel('Normalized FR Slope')
                    set(gca, 'ylim', [1.5*(min(animal_structure(id).FRslope(2,:))), 1.5*(max(animal_structure(id).FRslope(2,:)))])
                    hold on
                    set (gca, 'fontsize', 15, 'box', 'off')
                    hold on
                    yyaxis right
                    hold on
                    line([min(animal_structure(id).FRslope(2,:)), max(animal_structure(id).FRslope(2,:))], [.7, .7],'Color', [.1,.1,.1], 'LineStyle', '--', 'LineWidth', .3)
                    plot(x_md, animal_structure(id).kde, 'r-')
                    hold on
                    %plot (xaxis_MD_slope_days, kde2, 'c')

                    ylabel('Normalized Probability')
                    set (gca, 'ylim', [0,1], 'box', 'off')
                    hold off
                    set (gca, 'fontsize', 15)
                    set (gca, 'xlim', [0,7])
            %}
            
end

