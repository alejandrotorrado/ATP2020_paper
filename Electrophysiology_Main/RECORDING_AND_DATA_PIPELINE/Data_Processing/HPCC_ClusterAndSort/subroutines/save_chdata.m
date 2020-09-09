function save_chdata(newblock,thischan,saveDir,autoclustdir)

chdata.block = newblock;
chdata.channel = thischan;

fprintf('Saving data for channel %u...\n',thischan);
savet0 = tic;
save([saveDir filesep autoclustdir filesep 'channel_' num2str(thischan) '.mat'],...
    'chdata','-v7.3');

savet1 = toc(savet0);
fprintf('That took %.2f seconds.\n\n',savet1);


end