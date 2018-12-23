c = categorical({'16procs/4nodes', '16procs/2nodes' , '8procs/4nodes' , '8procs/2nodes' , '4procs/4nodes','4procs/2nodes','2procs/2nodes'});
c = reordercats(c ,{'16procs/4nodes', '16procs/2nodes' , '8procs/4nodes' , '8procs/2nodes', '4procs/4nodes','4procs/2nodes','2procs/2nodes'});
Times = [14-214131e-6 , 17.114123 , 13 - 0.116194 , 11-0.405144 , 16.393554 ,18.140432,21]; 
TimesKnn = [4.251333 ,  5 - 0.182228  , 9.189922 , 7 - 0.113627 , 19 - 0.457924, 23 - 0.030620,32-0.129771]; 

b = bar(c,Times,'FaceColor','flat');

figure(2)
b3 = bar(c,TimesKnn,'FaceColor','flat');

figure(3)
b2 = bar(c,[Times;TimesKnn]','stacked')
