function calculateBatchSpectras(all_names, startFreq, stopFreq, ...
    discretization, sAlpha, sBeta, ofrequencies, oAIP, oF)

frequencies = linspace(startFreq,stopFreq,(stopFreq-startFreq)/discretization +1);

    for i=1:length(all_names)
        filename = all_names(i);
        file = open(filename+string('.mat'));
        xval = file.xval;
        x_history = file.x_history;
        fobj = file.fobj;
        fres = file.fres;
        io = file.data_struct;
        clear file
        t = (io.tmax-io.tmin)*xval/100 +io.tmin;
        [K, M] = shellMatricesFromElements(io.Bke, io.Ske, io.Me, io.lm, t, io.YOUNG, io.RHO);
        C = sAlpha*M +sBeta*K;
        FF = oF{io.force_type};
        dr_dofs = io.dr_dofs;
        free_dofs = io.free_dofs;
        dr_values = io.dr_values;
        AIP = zeros(length(frequencies),1);
        parfor j=1:length(frequencies)
            omega = 2*pi*frequencies(j);
            Kd = K +1j*omega*C -omega*omega*M;
            u = SolveDirichletSystem(Kd, FF, dr_dofs, free_dofs, dr_values);
            AIP(j) = real(0.5*omega*omega*u'*C*u);
        end

    [sz1, sz2] = size(x_history);
    nsub = sqrt([sz1 sz1]);
    figure
        a2 = subplot(2,2,1);
        imagesc(reshape(xval,nsub))
        daspect([1 1 1])
        colormap(jet)
        colorbar
        a2.XAxis.Visible = 'off';
        a2.YAxis.Visible = 'off';
        title('Optimized Topology','FontWeight','bold')
        a3 = subplot(2,2,2);
        plot((fobj),'LineWidth',2)
        hold on
        plot(fres(:,1),'LineWidth',2)
        plot(fres(:,2),'LineWidth',2)
        grid on
        title('Objective Function and Restrictions History','FontWeight','bold')
        xlabel('Iteration #','FontWeight','bold')
        a1 = subplot(2,2,[3 4]);
        semilogy(ofrequencies,oAIP{io.force_type},'LineWidth',2)
        hold on
        semilogy(frequencies,AIP,'b','LineWidth',2)
        xline(io.freq,'--k',{'Excitation'},'FontWeight','bold','FontSize',10)
        xlim([startFreq stopFreq])
        legend('Initial topology', 'Optimized topology')
        title('Active Input Power Spectrum','FontWeight','bold')
        xlabel('Frequency [Hz]','FontWeight','bold')
        name = strcat(filename,'.png');
        saveas(gcf,name);
    end
end