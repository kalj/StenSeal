function compare_cartesian_general( benchmark, location) 

run(benchmark);

UpwindCartesian = {
    Second_order_Upwind_Cartesian , '2nd o. C';
    Fourth_order_Upwind_Cartesian , '4th o. C';
    Sixth_order_Upwind_Cartesian , '6th o. C';
};

UpwindGeneral = {
    Second_order_Upwind_General ,'2nd o. G';
    Fourth_order_Upwind_General,'4th o. G';
    Sixth_order_Upwind_General ,'6th o. G';
};


%% Orders
%% Matrix/Matrix free
%% Upwind/Compact
%% Cartesian/General


m = [UpwindCartesian{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, UpwindCartesian(:,2), UpwindCartesian(:,1), UpwindGeneral(:,2), UpwindGeneral(:,1));
title('Cartesian vs General, upwind','interpreter','latex');

savepng(fh, [location '\cartesian_vs_general_upwind'])

CompactCartesian = {
    Second_order_Compact_Cartesian, '2nd o. C';
    Fourth_order_Compact_Cartesian, '4th o. C';
    Sixth_order_Compact_Cartesian,  '6th o. C';
};


CompactGeneral = {
    Second_order_Compact_General, '2nd o. G';
    Fourth_order_Compact_General, '4th o. G';
    Sixth_order_Compact_General,  '6th o. G';
};



%% Orders
%% Matrix/Matrix free
%% Upwind/Compact
%% Cartesian/General


m = [UpwindCartesian{1}{:,1}];
[fh, ph, lh] = dofs_time_plot(m, CompactCartesian(:,2), CompactCartesian(:,1), CompactGeneral(:,2), CompactGeneral(:,1));
title('Cartesian vs General, Compact, General','interpreter','latex');

savepng(fh, [location '\cartesian_vs_general_compact'])

end
