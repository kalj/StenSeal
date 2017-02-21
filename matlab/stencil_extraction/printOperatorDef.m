function printOperatorDef()
    m = 25;
    nDig = 25;

    ops = {@op2nd, @op4th, @op6th};

    for i = 1:length(ops)
        op = ops{i}(m);
        h = op.h;
        c = op.c;
        Hv = subs(op.Hv, h, 1);
        e_l = subs(op.e_l, h, 1);
        d1_l = subs(op.d1_l, h, 1);
        D1_inner = subs(op.D1_inner, h, 1);
        D1_left_boundary = subs(op.D1_left_boundary, h, 1);
        D1_right_boundary = subs(op.D1_right_boundary, h, 1);
        D2_boundary = subs(op.D2_boundary, h, 1);
        D2 = subs(op.D2, h, 1);
        D2 = collect(D2, [c, h]);
        D2_boundary = collect(D2_boundary, [c, h]);

        % H
        fprintf('H: %s\n\n', char(Hv));

        % e_l, d_l
        fprintf(' e_l: %s\n', char(e_l));
        fprintf('d1_l: %s\n\n', char(d1_l));

        % D1
        fprintf('D1 inner:\n')
        n = length(D1_inner);
        cs = vectorToCellStencil(D1_inner, (n+1)/2, nDig);
        fprintf('%s\n\n', cellStencilToStencil(cs, 'd1_interior'));

        fprintf('D1 left boundary:\n')
        D1_left_boundary = trimMatrix(D1_left_boundary);
        cs = {};
        for k = 1:size(D1_left_boundary,1)
            cs(k,:) = vectorToCellStencil(D1_left_boundary(k,:), k, nDig);
        end
        fprintf('%s\n\n', cellStencilToStencilTensor(cs, 'd1_boundary'));

        fprintf('D1 right boundary:\n')
        D1_right_boundary = trimMatrix(D1_right_boundary);
        cs = {};
        d = size(D1_right_boundary,2) - size(D1_right_boundary,1);
        for k = 1:size(D1_right_boundary,1)
            cs(k,:) = vectorToCellStencil(D1_right_boundary(k,:), k+d, nDig);
        end

        fprintf('%s\n\n', cellStencilToStencilTensor(cs, 'd1_boundary_r'));

        % D2
        fprintf('D2 inner:\n')
        j = floor(m/2);
        D2_inner = D2(j,:);

        cs = symbolicStencilToCellStencilTensor(D2_inner, c, j, nDig);
        fprintf('%s\n\n', cellStencilToStencilTensor(cs, 'd2_interior'));


        fprintf('D2 boundary:\n')
        fprintf('%s\n\n', D2boundaryToString(D2_boundary, c, nDig));
    end
end


function str = cellStencilToStencil(c, name)
    n = length(c);
    stencil = strjoin(c, ' + ');
    str = sprintf('constexpr stenseal::Stencil<%d> %s(%s);', n, name, stencil);
end

function str = cellStencilToStencilTensor(c, name)
    [n,m] = size(c);

    typename = sprintf('constexpr stenseal::StencilTensor2D<%d,%d> %s', n,m,name);

    prefix = sprintf('%s(', typename);
    suffix = sprintf(');\n');
    spacer = repmat(' ', [1, length(prefix)]);

    % Calculate column width
    colWidth = {};
    for j = 1:m
        colWidth{j} = 0;
        for i = 1:n
            colWidth{j} = max(colWidth{j}, length(c{i,j}));
        end
    end

    % Pad all the strings
    for i = 1:n
        for j = 1:m
            c{i,j} = sprintf('%*s', colWidth{j}, c{i,j});
        end
    end

    % combine the stencils
    s = {};
    for i = 1:size(c,1)
        s{i} = strjoin(c(i,:), ' + ');
    end

    stencils = strjoin(s, sprintf(',\n%s', spacer));

    str = [prefix, stencils, suffix];
end

% convert the vector v to a stencil around index j
function c = vectorToCellStencil(v, j, n)
    c = cell(1, length(v));
    for i = 1:length(v)
        c{i} = sprintf('(%s)*sym[%d]', sym2string(v(i), n), i-j);
    end
end

% Convert a symbolic number to a string.
% Make sure the string representation is shorten than n chars.
function str = sym2string(val, n)
    p = n-5;
    assert(p > 15);

    str = char(vpa(val, p));

    % We could try to preserve the fractions but then we need to put .0 everywhere to aviod integer division.

    % str = char(val);

    % if len(str) > n
    %     p = n-5;
    %     assert(p > 15);

    %     str = char(vpa(val, p));
    % end
end

function C = symbolicStencilToCellStencilTensor(symSten, c, j, nDig)
    V = sym([]);
    for i = 1:length(symSten)
        V(i,:) = notStupidCoeffs(symSten(i),c);
    end

    [I,J] = find(V);

    irange = min(I):max(I);
    jrange = min(J):max(J);

    j = j - (min(I) - 1);
    V = V(irange, jrange);

    C = {};
    for i = 1:size(V,1)
        C(i,:) = vectorToCellStencil(V(i,:), j, nDig);
    end
end

function B = trimMatrix(A)
    [I,J] = find(A);

    irange = min(I):max(I);
    jrange = min(J):max(J);

    B = A(irange, jrange);
end


function str = D2boundaryToString(D2_boundary, c, nDig)
    D2_boundary = trimMatrix(D2_boundary);

    c_stencils = {};
    for i = 1:size(D2_boundary, 1)
        for j = 1:size(D2_boundary, 2)
            c_stencils{i,j} = notStupidCoeffs(D2_boundary(i,j), c);
        end
    end

    % Find c range
    imin = inf;
    imax = 0;

    for i = 1:size(D2_boundary, 1)
        for j = 1:size(D2_boundary, 2)
            I = find(c_stencils{i,j});
            if isempty(I)
                continue
            end
            imin = min(imin, min(I));
            imax = max(imax, max(I));
        end
    end

    % Shorten all the c stencils
    for i = 1:size(D2_boundary, 1)
        for j = 1:size(D2_boundary, 2)
            c_stencils{i,j} = c_stencils{i,j}(imin:imax);
        end
    end


    boundaryStencils = {};
    for i = 1:size(c_stencils,1)
        boundaryStencils{i} = {};
        for j = 1:size(c_stencils,2)
            boundaryStencils{i}(j,:) = vectorToCellStencil(c_stencils{i,j}, j, nDig);
        end
    end

    str = '';
    for i = 1:length(boundaryStencils)
        bdry = cellStencilToStencilTensor(boundaryStencils{i}, sprintf('bdry_row%d', i));
        str = [str sprintf('%s\n', bdry)];
    end


    k = length(boundaryStencils);
    [n,m] = size(boundaryStencils{1});
    tensor3dString = sprintf('constexpr stenseal::StencilTensor3D<%d,%d,%d> d2_boundary(', k, n, m);

    spacer = repmat(' ', [1, length(tensor3dString)]);

    contentCell = {};
    for i = 1:k
        contentCell{i} = sprintf('bdry_row%d', i);
    end

    content = strjoin(contentCell, sprintf(',\n%s',spacer));
    str = [str, sprintf('%s%s);\n', tensor3dString, content)];
end
