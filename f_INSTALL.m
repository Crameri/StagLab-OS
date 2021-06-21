
%%                                                     StagLab INSTALL 1.4

%                                                Fabio Crameri, 29.06.2019
function f_INSTALL

mFilesInUse         = dbstack('-completenames');
currentParfile      = mFilesInUse(end,1).file; %original file (i.e., parfile)
currentParfilePath  = currentParfile(1,1:max(strfind(currentParfile,filesep))); %path to original file

cd(currentParfilePath);
if exist(['.',filesep,'StagLab'],'dir')
    cd(['.',filesep,'StagLab']);            %go to the StagLab folder
elseif exist(['..',filesep,'StagLab'],'dir')
    cd(['..',filesep,'StagLab']);           %go to the StagLab folder
elseif exist(['..',filesep,'..',filesep,'StagLab'],'dir')
    cd(['..',filesep,'..',filesep,'StagLab']);        %go to the StagLab folder
else
    error('Could not find StagLab directory')
end

fAIo.task   = 'installStagLab';
[~,~] = f_AIo(fAIo,[],[],[],[],[]);
