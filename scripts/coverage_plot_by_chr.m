FONT_SIZE=20;

CHROMOSOMES={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"};
CHROMOSOME_LENGTHS=load('chr_lengths.txt');

figure(1); hold on;
for i=[1:length(CHROMOSOMES)]
    A=load(sprintf('ST001_counts_%s.csv',CHROMOSOMES{i}));
    h=plot(A(:,1),A(:,2)./CHROMOSOME_LENGTHS(i),'.-');
    if (i==16 || i==21)
        set(h,'LineWidth',2);
    endif
endfor
legend(CHROMOSOMES,'location','northeastoutside');
xlabel('Coverage (X)'); ylabel('Number of raw calls per base'); title('ST001 - PacBio');
axis square; grid on; set(gca,'fontsize',FONT_SIZE);
