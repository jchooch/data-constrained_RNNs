%% csv_writer.m
% This converts ZNN calcium activities from Matlab data format to CSV

bio_data = load('bio_data.mat'); % this loads the calcium activities (the biological data) as a structure
znn_acts = bio_data.fr;    % this gets the actual data matrix from the structure
csvwrite('znn_acts.csv', znn_acts); % this writes the data matrix of znn activities to a csv file

% The below code converts the stimuli metadata to csv format too (I should get the original stimuli metadata though...) 

[tLR,tLU,tLL,tLD,tRR,tRU,tRL,tRD] = deal(zeros(1,length(t)));   % Whit hard-coded all of this based on his reading of Matt's metadata file. Must check!
tLR(1184:1214)=1; tLR(1752:1782)=1; tLR(2774:2804)=1;
tLU(162:190)=1; tLU(502:532)=1; tLU(730:760)=1; tLU(1070:1100)=1; tLU(1638:1668)=1; tLU(2206:2236)=1; tLU(2548:2578)=1; tLU(2662:2690)=1; tLU(2888:2918)=1; tLU(3116:3146)=1;
tLL(48:78)=1;tLL(1980:2010)=1;tLL(2320:2350)=1;tLL(3002:3032)=1;
tLD(274:304)=1;tLD(842:872)=1;tLD(1412:1440)=1;tLD(2434:2464)=1;tLD(3230:3260)=1;tLD(3684:3714)=1;
tRR(1184:1214)=1;tRR(1752:1782)=1;tRR(2774:2804)=1;
tRU(502:532)=1;tRU(616:646)=1;tRU(956:986)=1;tRU(1070:1100)=1;tRU(1298:1328)=1;tRU(1412:1440)=1;tRU(1524:1554)=1;tRU(2548:2578)=1;tRU(2662:2690)=1;
tRL(48:78)=1;tRL(1980:2010)=1;tRL(2320:2350)=1;tRL(3002:3032)=1;
tRD(274:304)=1;tRD(730:760)=1;tRD(1638:1668)=1;tRD(1866:1896)=1;tRD(2092:2122)=1;tRD(2206:2236)=1;tRD(3116:3146)=1;tRD(3342:3372)=1;tRD(3570:3600)=1;tRD(3798:3828)=1;tRD(3912:3940)=1;

% Create an 8x4001 (stimuli types x timepoints) stimcourse metadata matrix
stimcourse = [tLR; tLU; tLL; tLD; tRR; tRU; tRL; tRD];  % Row order: LR,LU,LL,LD,RR,RU,RL,RD (left eye{right,up,left,down},right eye{right,up,left,down})

csvwrite('stimcourse.csv', stimcourse) % this writes the stimcourse metadata matrix to a csv file
    
%% scrap
%{
For k = 1:100
    plot(real(S.(sprintf('Y%d',k))))
    
Files=dir('*.mat');
N=length(Files);
Data=cell(N,1);
for k=1:N
Data{k}=load(Files(k).name);
end
%}
