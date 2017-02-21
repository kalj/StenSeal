data;

matrixFreeUpwind = {
    matrixFree_Second_order_Upwind, '2nd o. MF';
    matrixFree_Fourth_order_Upwind, '4th o. MF';
};

matrixUpwind = {
    matrix_Second_order_Upwind,'2nd o. M';
    matrix_Fourth_order_Upwind,'4th o. M';
};

matrixFreeCompact = {
    matrixFree_Second_order_Compact, '2nd o. MF';
    matrixFree_Fourth_order_Compact, '4th o. MF';
};

matrixCompact = {
    matrix_Second_order_Compact,'2nd o. M';
    matrix_Fourth_order_Compact,'4th o. M';
};


%% Orders
%% Matrix/Matrix free
%% Upwind/Compact
%% Cartesian/General


m = [matrixFreeCompact{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, matrixFreeCompact(:,2), matrixFreeCompact(:,1), matrixCompact(:,2), matrixCompact(:,1));
title('Matrix vs. MatrixFree, Compact','interpreter','latex');

savepng(fh, 'matrix_vs_matrixfree_compact')


m = [matrixFreeUpwind{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, matrixFreeUpwind(:,2), matrixFreeUpwind(:,1), matrixUpwind(:,2), matrixUpwind(:,1));
title('Matrix vs. MatrixFree, Upwind','interpreter','latex');

savepng(fh, 'matrix_vs_matrixfree_upwind')
