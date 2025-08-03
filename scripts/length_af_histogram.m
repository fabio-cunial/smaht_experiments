FONT_SIZE=16;

A=load('ST001_230_lengths.csv');
B=load('ST001_230_support_vaf.csv');


# SVLEN
figure(1); hold on;

bar([1:length(A)],A(:,[4,3]),'grouped');
set(gca, 'xtick', [1:length(A)]);
set(gca, 'xticklabel',{'100', '200', '300', '400', '500', '600', '700', '800', '900', '1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k'});
xlabel('SVLEN'); ylabel('Number of raw calls'); title('ST001 - PacBio 230x');
legend('DEL','INS','location','northeast'); grid on; set(gca,'fontsize',FONT_SIZE);


# Support and AF
figure(2); 

subplot(1,2,1); hold on;
hist(B(:,1),100);
xlabel('Number of supporting reads'); ylabel('Number of raw calls'); title('ST001 - PacBio 230x');
axis square; grid on; set(gca,'fontsize',FONT_SIZE);

subplot(1,2,2); hold on;
hist(B(:,2),50);
xlabel('VAF'); ylabel('Number of raw calls'); title('ST001 - PacBio 230x');
axis square; grid on; set(gca,'fontsize',FONT_SIZE);
