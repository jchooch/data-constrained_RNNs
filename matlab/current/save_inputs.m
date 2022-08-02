% save_inputs.m

T = 4001;
C = 8;

inputs = zeros(C, T);

inputs(1, 1184:1214) = 1;   % Right Up (REUM) V
inputs(1, 1752:1782) = 1;
inputs(1, 2774:2804) = 1;

inputs(2, 162:190) = 1;     % Right Right (RERM) V
inputs(2, 502:532) = 1;
inputs(2, 730:760) = 1;
inputs(2, 1070:1100) = 1;
inputs(2, 1638:1668) = 1;
inputs(2, 2206:2236) = 1;
inputs(2, 2548:2578) = 1;
inputs(2, 2662:2690) = 1;
inputs(2, 2888:2918) = 1;
inputs(2, 3116:3146) = 1;

inputs(3, 48:78) = 1;       % Right Down (REDM) V
inputs(3, 1980:2010) = 1;
inputs(3, 2320:2350) = 1;
inputs(3, 3002:3032) = 1;

inputs(4, 274:304) = 1;     % Right Left (RELM) V
inputs(4, 842:872) = 1;
inputs(4, 1412:1440) = 1;
inputs(4, 2434:2464) = 1;
inputs(4, 3230:3260) = 1;
inputs(4, 3684:3714) = 1;

inputs(5, 1184:1214) = 1;   % Left Up (LEUM) V
inputs(5, 1752:1782) = 1;
inputs(5, 2774:2804) = 1;
inputs(6, 502:532) = 1;     % Left Right (LERM) V
inputs(6, 616:646) = 1;
inputs(6, 956:986) = 1;
inputs(6, 1070:1100) = 1;
inputs(6, 1298:1328) = 1;
inputs(6, 1412:1440) = 1;
inputs(6, 1524:1554) = 1;
inputs(6, 2548:2578) = 1;
inputs(6, 2662:2690) = 1;
inputs(7, 48:78) = 1;       % Left Down (LEDM) V
inputs(7, 1980:2010) = 1;
inputs(7, 2320:2350) = 1;
inputs(7, 3002:3032) = 1;
inputs(8, 274:304)=1;       % Left Left (LELM) V
inputs(8, 730:760)=1;
inputs(8, 1638:1668)=1;
inputs(8, 1866:1896)=1;
inputs(8, 2092:2122)=1;
inputs(8, 2206:2236)=1;
inputs(8, 3116:3146)=1;
inputs(8, 3342:3372)=1;
inputs(8, 3570:3600)=1;
inputs(8, 3798:3828)=1;
inputs(8, 3912:3940)=1;

save('inputs.mat','inputs');

%tstim

tRstim = zeros(1,2000);
tRstim(:,592:607)=1;tRstim(:,876:891)=1;tRstim(:,1387:1402)=1;
tRstim(:,81:96)=2;tRstim(:,251:266)=2;tRstim(:,365:380)=2;tRstim(:,535:550)=2;tRstim(:,819:834)=2;tRstim(:,1103:1118)=2;tRstim(:,1273:1288)=2;tRstim(:,1330:1345)=2;tRstim(:,1443:1458)=2;tRstim(:,1557:1572)=2;
tRstim(:,24:39)=3;tRstim(:,990:1005)=3;tRstim(:,1160:1175)=3;tRstim(:,1501:1516)=3;
tRstim(:,137:152)=4;tRstim(:,421:436)=4;tRstim(:,706:721)=4;tRstim(:,1217:1232)=4;tRstim(:,1615:1630)=4;tRstim(:,1842:1857)=4;

tLstim = zeros(1,2000);
tLstim(:,592:607)=1;tLstim(:,876:891)=1;tLstim(:,1387:1402)=1;
tLstim(:,251:266)=2;tLstim(:,308:323)=2;tLstim(:,478:493)=2;tLstim(:,535:550)=2;tLstim(:,649:664)=2;tLstim(:,706:721)=2;tLstim(:,762:777)=2;tLstim(:,1273:1288)=2;tLstim(:,1330:1345)=2;
tLstim(:,24:39)=3;tLstim(:,990:1005)=3;tLstim(:,1160:1175)=3;tLstim(:,1501:1516)=3;
tLstim(:,137:152)=4;tLstim(:,365:380)=4;tLstim(:,819:834)=4;tLstim(:,933:948)=4;tLstim(:,1046:1061)=4;tLstim(:,1103:1118)=4;tLstim(:,1557:1572)=4;tLstim(:,1671:1686)=4;tLstim(:,1785:1800)=4;tLstim(:,1899:1914)=4;tLstim(:,1955:1970)=4;

left_stim_times = tLstim;
right_stim_times = tRstim;

save('left_stim_times.mat','left_stim_times');
save('right_stim_times.mat','right_stim_times');
