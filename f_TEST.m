
%%                                                        TEST StagLab 2.3

%                                                Fabio Crameri, 03.05.2020
clear
clf

StagLabDirectory    = which('SL_FieldPlot');
StagLabDirectory    = erase(StagLabDirectory,'SL_FieldPlot.m');
NumberTests         = 4;
TestNumber          = 0;
cd(fullfile(StagLabDirectory,'Examples','ExampleSetups'))
originalDirectory   = pwd;

%% START UP 
wb = waitbar(0,'Please wait...');
%wb = waitbar(0,'Please wait...','WindowStyle', 'modal');
wbPos       = get(wb,'Position');
screenSize  = get(groot,'Screensize');
set(wb,'Position',[screenSize(3)-wbPos(3),screenSize(4)-wbPos(4),wbPos(3),wbPos(4)]);
if ~exist(fullfile('..','ExampleFigures'),'dir')
    mkdir(fullfile('..','ExampleFigures'));
end


%% TEST 1
TestNumber          = TestNumber+1;
try
    [TEST] = SetupStagYYndSubduction;
catch me
    TEST.successful = false;
end
if TEST.successful
    string          = 'successful';
else
    string          = 'failed';
    disp(' '); disp(me.message)
end
disp(' ')
disp(['Test ',num2str(TestNumber),' (StagYY,binary,non-dimensional,Subduction)'])
disp(['  ',char(10551),' ',string])
cd(originalDirectory)
TestConclusion(TestNumber,1) 	= TEST.successful;
waitbar(TestNumber/NumberTests,wb)


%% TEST 2
TestNumber          = TestNumber+1;
try
    [TEST] = SetupStagYYdimVenus;
catch me
    TEST.successful = false;
end
if TEST.successful
    string          = 'successful';
else
    string          = 'failed';
    disp(' '); disp(me.message)
end
disp(' ')
disp(['Test ',num2str(TestNumber),' (StagYY,binary,dimensional,Venus)'])
disp(['  ',char(10551),' ',string])
cd(originalDirectory)
TestConclusion(TestNumber,1) 	= TEST.successful;
waitbar(TestNumber/NumberTests,wb)


%% TEST 3
TestNumber          = TestNumber+1;
try
    [TEST] = SetupStagYYtimedat;
catch me
    TEST.successful = false;
end
if TEST.successful
    string          = 'successful';
else
    string          = 'failed';
    disp(' '); disp(me.message)
end
disp(' ')
disp(['Test ',num2str(TestNumber),' (StagYY,binary,dimensional,timedat)'])
disp(['  ',char(10551),' ',string])
cd(originalDirectory)
TestConclusion(TestNumber,1) 	= TEST.successful;
waitbar(TestNumber/NumberTests,wb)


%% TEST 4
TestNumber          = TestNumber+1;
try
    [TEST] = SetupStagYYrprof;
catch me
    TEST.successful = false;
end
if TEST.successful
    string          = 'successful';
else
    string          = 'failed';
    disp(' '); disp(me.message)
end
disp(' ')
disp(['Test ',num2str(TestNumber),' (StagYY,binary,dimensional,rprof)'])
disp(['  ',char(10551),' ',string])
cd(originalDirectory)
TestConclusion(TestNumber,1) 	= TEST.successful;
waitbar(TestNumber/NumberTests,wb)


%% FINISHING UP
close all
disp(' '); disp(' ');
if min(TestConclusion==0)
    disp([char(9888),' Unfortunately, some problems occurred.',char(9888)])
else
    disp([char(9745),' All Tests Successful.'])
end
close(wb)

