% close all
% figure(1)
    hold on;
    triplot(Elements, Nodes(:,1),Nodes(:,2));
    xlabel('X','fontsize',14)
    ylabel('Y','fontsize',14)
    title('Mesh','fontsize',14)
    fh = figure(1);
    hold off
%    set(fh, 'color', 'white');
    
%    temp = Edges(Boundary, 1:2);
%    scatter(Nodes(temp(:), 1), Nodes(temp(:), 2), 18, 'k', 'filled');
    
%    [dim, dmp] = size(Global_DoF_Addrs);
%    
%    textCell = arrayfun(@(n) sprintf('(GDoF, %d)',n),[1:dim]','un',0);
%    
%    scatter(Global_DoF_Addrs(:,1), Global_DoF_Addrs(:,2), 12, 'r');
    % Add textCell
%for ii = 1:dim 
%    text(Global_DoF_Addrs(ii,1)+.02, Global_DoF_Addrs(ii,2)+.02,textCell{ii,1},'FontSize',12) 
%end
%    hold off;
%    
%    ax = gca;
%    set (ax, "visible", "off");
%    print (ax, "Sample_Global_Mesh.pdf", "-dpdflatexstandalone");
%    set (ax, "visible", "on");