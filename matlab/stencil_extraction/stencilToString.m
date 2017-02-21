% Takes a stencil for point j applied to the symbolic vector c and converts it
% To stenseal form.
function str = stencilToString(j, stencil, c)
    str = '';
    co = notStupidCoeffs(stencil, c);
    I = find(co);
    for i = 1:length(I)-1
        str = [str, sprintf('%s*sym[%d] + ', char(co(I(i))), I(i)-j)];
    end
    str = [str, sprintf('%s*sym[%d],',char(co(I(end))), I(end)-j)];
end