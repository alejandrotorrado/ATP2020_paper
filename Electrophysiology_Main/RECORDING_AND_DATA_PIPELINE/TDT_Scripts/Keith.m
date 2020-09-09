Tank = 'C:\TDT\OpenEx\Tanks\DemoTank2';                     %Tank Name and Location
Block = 'Block-13';                                         %Block Name

TT = actxcontrol('TTank.X');
TT.ConnectServer('Local','Me');
TT.OpenTank(Tank, 'R');
TT.SelectBlock(Block);

notes = TT.CurBlockNotes;                                   %Block Notes
duration =  TT.CurBlockStopTime-TT.CurBlockStartTime;       %Recording Duration in Seconds

TT.CloseTank;
TT.ReleaseServer;