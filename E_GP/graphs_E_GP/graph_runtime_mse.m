% graph_runtime_mse draw the comparison of two methods: srch, dpstrf(cholesky with diagonal pivoting),

% x-axis is target rank
% y-axis is the mse(mean squared error)

rank = [ 72
         103
         135
         186
         326
         404
         454
         524
         553
         604
         743
         819
         994];
srch = [  
  0.24102599999999996 
  0.26968899999999985     
  0.28319700000000014  
  0.32802199999999981 
  0.38528399999999996     
  0.45350299999999999     
  0.45592699999999997     
  0.54625599999999985
  0.56035299999999988
  0.57138999999999984          
  0.78092099999999998     
  0.85298800000000030     
   1.0389579999999998 ];
dpstrf = [ 
  0.12650899999999998
  0.13297900000000018 
  0.23424999999999985     
  0.24778100000000025     
  0.56537999999999977     
  0.68118699999999999     
  0.78366999999999987     
  0.81891200000000008     
  0.89157600000000059     
  0.98886199999999969     
   1.1780910000000002     
   1.1778150000000003     
   1.4280819999999999  ];

srch_prediction = [
24.366557781395610
17.721917411985668
16.984182662512424
16.747641519512559
16.485487565532534
16.466088257330764
16.403508917722601
16.382124579591881 
16.357301789893828
16.380405785670991
16.335056054615301
16.318311081414084
16.319736439524302];
dpstrf_prediction = [
   22.676194432342029
   18.265367130541701
   16.870564263539524     
   16.670861323955933     
   16.508985154562254     
   16.421186433825817     
   16.395671244236510     
   16.389995090458832     
   16.399219538329938     
   16.393076031567372     
   16.385856468555055     
   16.369522207663209        
   16.336230437288606 ]; 

figure(1)
plot(rank, dpstrf, '--*k', 'LineWidth',1.5);
hold on
plot(rank, srch, '--or', 'LineWidth',1.5);
hold off
leg = legend('DPSTRF', 'SRCH', 'Location', 'northwest');
set(leg, 'FontSize', 20)
xlabel('Approximate rank', 'FontSize', 20);
ylabel('Run time (sec)', 'FontSize', 20);
title('Run time vs. Approximate rank', 'FontSize', 20)

figure(2)
plot(rank, dpstrf_prediction, '--*k', 'LineWidth',1.5);
hold on
plot(rank, srch_prediction, '--or', 'LineWidth',1.5);
hold off
leg=legend('DPSTRF', 'SRCH', 'Location', 'northeast');
set(leg,'FontSize',20)
xlabel('Approximate rank','FontSize',20);
ylabel('Mean squared prediction error','Fontsize',20);
title('MSE vs. Approximation rank','Fontsize',20)