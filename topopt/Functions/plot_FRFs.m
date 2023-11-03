function plotFRFs(parameters,frequencies,FRF_init,FRF_final)
switch parameters.objective_function
     case {"v2_rms", "v2_db"}
         idx = 1;
         ystr = 'Quadratic Mean Velocity [m/s]';
    case "AIP"
        idx = 2;
        ystr = 'Active Input Power [W]';
    case "dynamic compliance"
        idx = 3;
        ystr = 'Dynamic Compliance [m/N]';
    otherwise
        ystr = 'Eigenvalue Gap';

end
    

loglog(frequencies,FRF_init(:,idx),'b','LineWidth',2)
hold on
loglog(frequencies,FRF_final(:,idx),'g','LineWidth',2)
grid on
xline(parameters.freq,'--k',{'Excitation Frequency'},'FontWeight','bold','FontSize',18)
legend('Initial Design','Optimized Design')
set(gca,'FontSize',20)
xlabel('Frequency [Hz]','FontSize',24,'FontWeight','bold')
ylabel(ystr,'FontSize',24,'FontWeight','bold')
title('Initial vs Optimized Design','FontWeight','bold','FontSize',24)
end