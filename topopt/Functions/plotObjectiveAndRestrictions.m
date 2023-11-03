function plotObjectiveAndRestrictions(fobj,fres,ystr)    
hold on
    subplot(1,2,1)
    semilogy(fobj,'LineWidth',2)
    grid on
    set(gca,'FontSize',20)
    xlabel('Iteration','FontSize',24,'FontWeight','bold')
    title('Objective Function History','FontWeight','bold','FontSize',24)
    ylabel(ystr,'FontSize',24,'FontWeight','bold')
    subplot(1,2,2)
    hold on
    for i=1:size(fres,2)
    semilogy(fres(:,i),'LineWidth',2)
    end
    grid on
    set(gca,'FontSize',20)
    legend('Restriction 1','Restriction 2')
    xlabel('Iteration','FontSize',24,'FontWeight','bold')
    title('Restriction Function History','FontWeight','bold','FontSize',24)
    ylabel('Mass Restriction','FontSize',24,'FontWeight','bold')
end
