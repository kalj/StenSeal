function compare_matrix_matrix_free( benchmark, location) 

run(benchmark);

matrixFreeUpwind = {
    Second_order_Upwind_Cartesian , '2nd o. MF';
    Fourth_order_Upwind_Cartesian , '4th o. MF';
    Sixth_order_Upwind_Cartesian , '6th o. MF';
};

matrixUpwind = {
    Second_order_Upwind_Cartesian_Matrix ,'2nd o. M';
    Fourth_order_Upwind_Cartesian_Matrix ,'4th o. M';
    Sixth_order_Upwind_Cartesian_Matrix ,'6th o. M';
};

matrixFreeCompact = {
    Second_order_Compact_Cartesian, '2nd o. MF';
    Fourth_order_Compact_Cartesian, '4th o. MF';
    Sixth_order_Compact_Cartesian,  '6th o. MF';
};

matrixCompact = {
    Second_order_Compact_Cartesian_Matrix,'2nd o. M';
    Fourth_order_Compact_Cartesian_Matrix,'4th o. M';
    Sixth_order_Compact_Cartesian_Matrix, '6th o. M'
};


%% Orders
%% Matrix/Matrix free
%% Upwind/Compact
%% Cartesian/General


m = [matrixFreeCompact{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, matrixFreeCompact(:,2), matrixFreeCompact(:,1), matrixCompact(:,2), matrixCompact(:,1));
title('Matrix vs. MatrixFree, Compact','interpreter','latex');

savepng(fh, [location '\matrix_vs_matrixfree_compact'])


m = [matrixFreeUpwind{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, matrixFreeUpwind(:,2), matrixFreeUpwind(:,1), matrixUpwind(:,2), matrixUpwind(:,1));
title('Matrix vs. MatrixFree, Upwind','interpreter','latex');

savepng(fh, [location '\matrix_vs_matrixfree_upwind'])


matrixFreeUpwind = {
    Second_order_Upwind_General , '2nd o. MF';
    Fourth_order_Upwind_General , '4th o. MF';
    Sixth_order_Upwind_General , '6th o. MF';
};

matrixUpwind = {
    Second_order_Upwind_General_Matrix ,'2nd o. M';
    Fourth_order_Upwind_General_Matrix ,'4th o. M';
    Sixth_order_Upwind_General_Matrix ,'6th o. M';
};

matrixFreeCompact = {
    Second_order_Compact_General, '2nd o. MF';
    Fourth_order_Compact_General, '4th o. MF';
    Sixth_order_Compact_General,  '6th o. MF';
};

matrixCompact = {
    Second_order_Compact_General_Matrix,'2nd o. M';
    Fourth_order_Compact_General_Matrix,'4th o. M';
    Sixth_order_Compact_General_Matrix, '6th o. M'
};


%% Orders
%% Matrix/Matrix free
%% Upwind/Compact
%% Cartesian/General


m = [matrixFreeCompact{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, matrixFreeCompact(:,2), matrixFreeCompact(:,1), matrixCompact(:,2), matrixCompact(:,1));
title('Matrix vs. MatrixFree, Compact, General','interpreter','latex');

savepng(fh, [location '\matrix_vs_matrixfree_compact_general'])


m = [matrixFreeUpwind{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, matrixFreeUpwind(:,2), matrixFreeUpwind(:,1), matrixUpwind(:,2), matrixUpwind(:,1));
title('Matrix vs. MatrixFree, Upwind, General','interpreter','latex');

savepng(fh, [ location '\matrix_vs_matrixfree_upwind_general'])
end
