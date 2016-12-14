% relative error comparison
index = 1:10;
dpstrf_s_20 = [0.4649, 0.4615, 0.5324, 0.6115, 0.6215, 0.6290, 0.6025, 0.6953, 0.7478, 0.7284];
srch_s_20 = [0.2963E-01,0.3684E-01,0.7554E-01,0.1068,0.1150,0.9812E-01,0.1606,0.1503,0.2180,0.1985] ;
dpstrf_s_40 = [0.1366, 0.1879, 0.3109, 0.2366, 0.2450, 0.2782, 0.2479, 0.3836, 0.3994, 0.4426];
srch_s_40 = [0.1050E-01,0.1104E-01,0.2125E-01,0.2313E-01,0.2335E-01,0.3148E-01,0.5872E-01,0.4118E-01,0.5734E-01,0.6782E-01];
dpstrf_s_60 = [0.9590E-01, 0.1444, 0.1598, 0.1649, 0.1515, 0.1735, 0.1985, 0.2490, 0.2327, 0.2361];
srch_s_60 = [0.5011E-02,0.5974E-02,0.7808E-02,0.1012E-01,0.1202E-01,0.1643E-01,0.3207E-01,0.1805E-01,0.3216E-01,0.3522E-01];
          
plot(index, dpstrf_s_20, '--+k', 'LineWidth',1.5);
hold on
plot(index, dpstrf_s_40, '--*c', 'LineWidth',1.5);
hold on
plot(index, dpstrf_s_60, '--sg', 'LineWidth',1.5);
hold on
plot(index, srch_s_20, '--+m', 'LineWidth',1.5);
hold on
plot(index, srch_s_40, '--*r', 'LineWidth',1.5);
hold on
plot(index, srch_s_60, '--sb', 'LineWidth',1.5);
hold off
leg=legend('DPSTRF,k=20','DPSTRF,k=40','DPSTRF,k=60','SRCH,k=20','SRCH,k=40','SRCH,k=60', 'Location', 'northwest');
set(leg,'FontSize',10);
xlabel('Top singular value index','FontSize', 20);
ylabel('Relative error','FontSize', 20);
title('Singular value approximation relative error','FontSize', 20)