vecCoords=readtable('/N/u/dnbulloc/Carbonate/Downloads/proj-5a9ee80853bd38003c02afa7/sub-100206/dt-neuro-dwi.tag-dtiinit.id-5c4f56e87f21220052ddc02b/dwi.csv')
vecArray = table2array( vecCoords )

figpath='/N/u/dnbulloc/Carbonate/Downloads/proj-5a9ee80853bd38003c02afa7/sub-100206/figs/spheres/'

for iAngles=1:length(vecArray)
    figure
    hold on
    quiver3(vecArray(1,:),vecArray(2,:),vecArray(3,:),vecArray(1,:)*2,vecArray(2,:)*2,vecArray(3,:)*2,0,'LineWidth',2)
    quiver3(vecArray(1,iAngles),vecArray(2,iAngles),vecArray(3,iAngles),vecArray(1,iAngles)*2,vecArray(2,iAngles)*2,vecArray(3,iAngles)*2,0,'LineWidth',7,'Color','r')
    view(45,45)
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    set(gca,'ztick',[],'zticklabel',[])
    axis equal
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];

    saveas(fig,strcat(figpath,'sphere_',num2str(iAngles),'.png'))
    close all
end


dwiTest=niftiRead('/N/u/dnbulloc/Carbonate/Downloads/proj-5a9ee80853bd38003c02afa7/sub-100206/dt-neuro-dwi.tag-dtiinit.id-5c4f56e87f21220052ddc02b/dwi.nii.gz')