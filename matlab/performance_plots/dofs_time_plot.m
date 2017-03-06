function [fh, ph, lh] = dofs_time_plot(m, group1legends, group1, group2legends, group2)
    fh = figure();

    ph = {};
    scale = 1e6;
    hold on
    for i = 1:length(group1)
        ph{1,i} = loglog(m, [group1{i}{:,2}]./m*scale);
        ph{1,i}.LineWidth = 2.0;
        ph{1,i}.Color = Color.colors{i};
    end

    for i = 1:length(group2)
        ph{1,i} = loglog(m, [group2{i}{:,2}]./m*scale);
        ph{1,i}.LineWidth = 2.0;
        ph{1,i}.LineStyle = '--';
        ph{1,i}.Color = Color.colors{i};
    end

    hold off

    lh = legend([group1legends; group2legends]);
    lh.Interpreter = 'latex';
    lh.Location = 'NorthEast';

    xlabel('N','interpreter','latex');
    ylabel('Time(ms)/(N 1e-6)', 'interpreter', 'latex')

    ah = gca;
    ah.TickLabelInterpreter = 'latex';

    setFontSize(fh);
end

