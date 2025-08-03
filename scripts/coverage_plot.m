FONT_SIZE=16;

B=load('SMHT001_counts.csv');
A=load('ST001_counts.csv');


figure(1); 

subplot(1,2,1); hold on;

plot(A(:,1),A(:,2),'.-k');
plot(A(:,1),A(:,3),'.-r');
plot(A(:,1),A(:,4),'.-b');
plot(A(:,1),A(:,5),'.-g');

plot(B(:,1),B(:,2),'ok');
plot(B(:,1),B(:,3),'or');
plot(B(:,1),B(:,4),'ob');
plot(B(:,1),B(:,5),'og');

xlabel('Coverage (X)'); ylabel('Number of raw calls'); title('ST001 - PacBio');
legend('All','INS','DEL','BND','location','northwest'); axis square; grid on; set(gca,'fontsize',FONT_SIZE);


subplot(1,2,2); hold on;

plot(A(:,1),A(:,6),'.-r');
plot(A(:,1),A(:,7),'.-k');

xlabel('Coverage (X)'); ylabel('Number of raw calls'); title('ST001 - PacBio');
legend('Inside TRs','Outside TRs','location','northwest'); axis square; grid on; set(gca,'fontsize',FONT_SIZE);
