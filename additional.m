% if plt_fig
%     figure, hold on, axis equal, xlim([0 1]), ylim([0 1]), zlim([0,1])
%     for i=1:size(T_Dspd,1)
%         plot3(ls(1,:,i),ls(2,:,i),m(1,:,i),'k')        
%         scatter3(ls(1,:,i),ls(2,:,i),ls(1,:,i)+ls(2,:,i),[],RGB(:,:,i)','v','filled')
%     end
%     xlabel('l'),ylabel('s'),zlabel('m');
%     
%     view(188,46)
%     
%     if plt_locus
%         MB_locus=LMSToMacBoyn(T_cones_sp);
%         %plot(MB_locus(1,:),MB_locus(2,:))
%         fill([MB_locus(1,5:65),MB_locus(1,5)],[MB_locus(2,5:65),MB_locus(2,5)],'k','LineStyle','none','FaceAlpha','0.1')
%     end
% end