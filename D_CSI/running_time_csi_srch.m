index = [50, 100, 150, 200];
data = zeros(4,4);

n = 3000; k = 50;
list1 = [6.921763871000000   7.659856471000000   5.290280742000000   4.162998550000000];
list2 = [7.003231807000000   8.405216894000000   5.929388669000000   4.134775478000000];
list3 = [7.014870083000000   8.334239997999999   5.470091516000000   4.521083742000000];
list4 = [7.031165561000000   8.159847518999999   5.458830141000000   4.478825381000000];
average = (list1 + list2 + list3 + list4) / 4.0;
data(1,:) = average;

n = 3000; k = 100;
list1 = [10.447410817000000  11.640220650000000   7.600420387000000   4.633115760000000];
list2 = [10.328210316000000  11.859193091000000   7.953464640000000   4.643220406000000];
list3 = [10.157720067000000  11.553616735000000   7.572077616000000   4.630005609000000];
list4 = [10.156614770999999  11.770634037000001   7.788106763000000   4.608770777000000];
average = (list1 + list2 + list3 + list4) / 4.0;
data(2,:) = average;

n = 3000; k = 150;
list1 = [13.839354548999999  15.656113154000000  10.244091481000000   5.234370951000000];
list2 = [13.981324970999999  15.846619006999999  10.380492428000000   5.430544162000000];
list3 = [14.545441278000000  16.156823502999998  10.637386738000000   5.201245797000000];
list4 = [13.822357215000000  15.671574644000000  10.308282223000001   5.185770061000000];
average = (list1 + list2 + list3 + list4) / 4.0;
data(3,:) = average;

n = 3000; k = 200;
list1 = [17.322935188999999  19.823429833999999  13.258310580000000   5.773350270000000];
list2 = [17.627835106999999  21.964475311000001  14.495564679999999   5.674940138000000];
list3 = [17.356996554999998  19.716839985000000  13.110072666000001   5.759203196000000];
list4 = [17.170043124999999  19.619617152000000  13.091188282999999   5.791090544000000];
average = (list1 + list2 + list3 + list4) / 4.0;
data(4,:) = average;

plot(index, data(:,2), '--+b', 'LineWidth',1.5);
hold on
plot(index, data(:,1), '--*g', 'LineWidth',1.5);
hold on
plot(index, data(:,3), '--sr', 'LineWidth',1.5);
hold on
plot(index, data(:,4), '--om', 'LineWidth',1.5);
hold off
leg=legend('CSI80','CSI40','Cholesky','SRCH','Location', 'northwest');
set(leg,'FontSize',20);
xlabel('Approximate rank','FontSize', 20);
ylabel('Run time (sec)','FontSize', 20);
title('Run time comparison on MNIST','FontSize', 20)