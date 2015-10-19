%Convert trc file to input file for vbSPT
clear all
%% Initialize fields (default vbSPT)
CylinderL = 40; % nm (length of cylindrical part only)
Radius = 20;  % nm    (spherical end caps, and cylinder radius)
% initiate options
timestep = 0; % [s]
stepSize = 5; %[nm]
locAccuracy = 0; %[nm]
transMat = 0; % [/timestep]
transRate = 0; % [/s]
occProb = 0;
Dapp = 0;
%trajLengths = 0;
runs = 1;
do_steadystate = false;
do_parallel = false;
do_single = false;
do_transRate = false;
do_transMat = false;
%Load datafile for conversion
[filename,filepath] = uigetfile('*.trc','Select trc file');
outputfile = sprintf('%s_converted.mat',filename);
output = fullfile(filepath, outputfile);
csvdata = csvread(filename);
numTraj = max(csvdata(:,1));
Traj = cell(1,numTraj); 
trajLengths = zeros(1,numTraj);
tnum = 1;
for i=1:length(csvdata)
    track = csvdata(i,1);
    frame = csvdata(i,2);
    x = csvdata(i,3);
    y = csvdata(i,4);
    if (i > 1 && track == cachetrack)    
        d = sqrt((x - cachex)^2 + (y-cachey)^2);
        timestep = frame - cacheframe;
    else
       d = 0;
       timestep = 1;
    end
    %cache
    cachex = x;
    cachey = y;
    cacheframe = frame;
    cachetrack = track;
    Traj{track} = cat(1,Traj{track},[x y round(timestep)]);
    %Traj{track}
    trajLengths(track) = trajLengths(track) + d;
    %Trajlengths(track)
end
%% List the parameters
X.trajLengths=trajLengths;
X.runs=runs;
X.do_steadystate=do_steadystate;
X.do_parallel=do_parallel;
X.do_single=do_single;
X.CylinderL=CylinderL;
X.Radius=Radius;
X.timestep=timestep;
X.stepSize=stepSize;
X.locAccuracy=locAccuracy;
X.finalTraj=Traj;
X.numTraj=numTraj; % = length(trajLengths)
X.avTrajLength=mean(trajLengths);
X.shortestTraj=min(trajLengths);
X.longestTraj=max(trajLengths);
X.Dapp=Dapp;
X.occProb=occProb;

save(output, '-struct','X');    
