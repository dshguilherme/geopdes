function calculateBatchSpectras(all_names, startFreq, stopFreq, ...
    discretization, sAlpha, sBeta)

frequencies = linspace(startFreq,stopFreq,(stopFreq-startFreq)/discretization +1);

    for i=1:length(all_names)
        filename = all_names(i);
        file = open(filename+string('.mat'));
        xval = file.xval;
        x_history = file.x_history;
        fobj = file.fobj;
        fres = file.fres;
        clear file
        io = data_struct;
        clear data_struct
        t = (io.tmax-io.tmin)*xval/100 +io.tmin;
        [K, M] = shellMatricesFromElements(io.Bke, io.Ske, io.Me, io.lm, t, io.YOUNG, io.RHO);
        C = sAlpha*K +sBeta*M;
        AIP = zeros(length(frequencies),1);
        parfor j=1:length(frequencies)
            omega = 2*pi*frequencies(j);
            Kd = K +1j*omega*C -omega*omega*M;
            u = SolveDirichletSystem(Kd, F, dr_dofs, free_dofs, dr_values);
            AIP(j) = real(0.5*omega*omega*u'*C*u);
        end

    [sz1, sz2] = size(x_history);
    nsub = sqrt([sz1 sz1]);
    figure
        a2 = subplot(2,1,1);
        imagesc(reshape(xval,nsub))
        daspect([1 1 1])
        colormap(jet)
        colorbar
        a2.XAxis.Visible = 'off';
        a2.YAxis.Visible = 'off';
        a1 = subplot(2,1,2);
        semilogy(ofrequencies,oAIP,'LineWidth',2)
        hold on
        semilogy(frequencies,AIP,'b','LineWidth',2)
        xline(io.freq,'--k',{'Excitation Frequency'},'FontWeight','bold','FontSize',18)
        xlim(frequencies)
        name = strcat(filename,'.png');
        saveas(gcf,name);
    end
end