figure(1);
plot(IRF2(end,:)-IRF0(end,:));
hold on;
plot(IRFLambda1(end,:)-IRF3(end,:));

figure(1);
plot(IRF2(end-1,:)-IRF0(end-1,:));
hold on;
plot(IRFLambda1(end-1,:)-IRF3(end-1,:));