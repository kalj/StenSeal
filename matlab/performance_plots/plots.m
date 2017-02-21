data;

generalGeomCases = {
    matrixFree_Second_order_Upwind, '2nd o. Upwind, MF';
    matrixFree_Fourth_order_Upwind, '4th o. Upwind, MF';
    matrixFree_Second_order_Compact, '2nd o. Compact, MF';
};

cartesianGeomCases = {
    matrixFree_Second_order_Upwind, '2nd o. Upwind, cartesian';
    matrixFree_Fourth_order_Upwind, '4th o. Upwind, cartesian';
    matrixFree_Second_order_Compact, '2nd o. Compact, cartesian';
};

%% Orders
%% Matrix/Matrix free
%% Upwind/Compact
%% Cartesian/General


fh = figure();
% ah = axes();

m = [Second_order_Compact{:,1}];

colors = {
    Color.blue,
    Color.red,
    Color.yellow,
    Color.purple,
    Color.green,
};

hold on
for i = 1:length(generalGeomCases)
    ph = loglog(m, [generalGeomCases{i}{:,2}]./m);
    ph.LineWidth = 2.0;
    ph.Color = colors{i};
end

for i = 1:length(cartesianGeomCases)
    ph = loglog(m, [cartesianGeomCases{i}{:,2}./m]);
    ph.LineWidth = 2.0;
    ph.LineStyle = '--';
    ph.Color = colors{i};
end
hold off

lh = legend([generalGeomCases(:,2); cartesianGeomCases(:,2)]);
lh.Interpreter = 'latex';
lh.Location = 'NorthWest';

title('Cartesian vs Curved geometry','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Time(ms)', 'interpreter', 'latex')

setFontSize(fh);

