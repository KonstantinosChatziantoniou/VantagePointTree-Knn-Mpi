%% data
c = categorical({'16procs/4nodes', '16procs/2nodes' , '8procs/4nodes' , '8procs/2nodes' , '8procs1node', '4procs/4nodes','4procs/2nodes',' 4procs1node' ,'2procs/2nodes','2procs1node' ,'1proc'});
c = reordercats(c ,{'16procs/4nodes', '16procs/2nodes' , '8procs/4nodes' , '8procs/2nodes' , '8procs1node', '4procs/4nodes','4procs/2nodes',' 4procs1node' ,'2procs/2nodes','2procs1node','1proc'});
VPtime = zeros(11,1);
Knntime = zeros(11,1);


%using fastest time
VPtime(1) = min([1.755594 , 1.755594 , 2.037593, 2.702081 ]);
Knntime(1) = min([1775 , 1291 , 1671 , 1775 ]);
VPtime(2) = min([1-0.065764, 3.641971 , 5-0.011990 , 1-0.065764 ]);
Knntime(2) = min([298 , 298 , 1638]);
VPtime(3) = min([1.860414, 2-0.052127 , 2 - 0.089894  ]);
Knntime(3) = min([1100 , 1104 , 981]);
VPtime(4) = min([2-0.376418 , 4.193004 , 2.292043 ]);
Knntime(4) = min([671 , 830 , 1042]);
VPtime(5) = min([4 - 0.644672 , 5 - 0.156545 , 2 - 0.478380]);
Knntime(5) = min([994 , 849 , 1106]);
VPtime(6) = min([4 - 0.490601, 3.5735 , 3.5533 ]);
Knntime(6) = min([667 , 710 , 553 ]);
VPtime(7) = min([3.544716 , 3.245795]);
Knntime(7) = min([730 , 482]);
VPtime(8) = min([5 - 0.4929, 3 - 0.365 , 3.4031]);
Knntime(8) = min([485 , 659 , 420]);
VPtime(9) = min([6.716 , 7 - 0.733 ]);
Knntime(9) = min([383 , 371]);
VPtime(10) = min([7 - 0.651771 ]);
Knntime(10) = min([347 ]);
VPtime(11) = [13.146];
Knntime(11) = [37-0.049];

meanVPtime = zeros(11,1);
meanKnntime = zeros(11,1);

meanVPtime(1) = mean([1.755594 , 1.755594 , 2.037593, 2.702081 ]);
meanKnntime(1) = mean([1775 , 1291 , 1671 , 1775 ]);
meanVPtime(2) = mean([1-0.065764, 3.641971 , 5-0.011990 , 1-0.065764 ]);
meanKnntime(2) = mean([298 , 298 , 1638]);
meanVPtime(3) = mean([1.860414, 2-0.052127 , 2 - 0.089894  ]);
meanKnntime(3) = mean([1100 , 1104 , 981]);
meanVPtime(4) = mean([2-0.376418 , 4.193004 , 2.292043 ]);
meanKnntime(4) = mean([671 , 830 , 1042]);
meanVPtime(5) = mean([4 - 0.644672 , 5 - 0.156545 , 2 - 0.478380]);
meanKnntime(5) = mean([994 , 849 , 1106]);
meanVPtime(6) = mean([4 - 0.490601, 3.5735 , 3.5533 ]);
meanKnntime(6) = mean([667 , 710 , 553 ]);
meanVPtime(7) = mean([3.544716 , 3.245795]);
meanKnntime(7) = mean([730 , 482]);
meanVPtime(8) = mean([5 - 0.4929, 3 - 0.365 , 3.4031]);
meanKnntime(8) = mean([485 , 659 , 420]);
meanVPtime(9) = mean([6.716 , 7 - 0.733 ]);
meanKnntime(9) = mean([383 , 371]);
meanVPtime(10) = 7 - 0.651771;
meanKnntime(10) = 347 ;
meanVPtime(11) = 13.146;
meanKnntime(11) = 37-0.049;

meanVPtime = meanVPtime - VPtime;
meanKnntime = meanKnntime - Knntime;

h = bar(c,[VPtime,meanVPtime],'stacked');
set(h,{'FaceColor'},{'c';'b'});
figure(2)
h2 = bar(c,[Knntime,meanKnntime],'stacked');
set(h2,{'FaceColor'},{'c';'b'});

