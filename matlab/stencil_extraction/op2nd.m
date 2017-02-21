function r = op2nd(m)
    h = sym('h', 'real');

    Hv = ones(m,1);
    Hv(1) = 1/2;
    Hv(m) = 1/2;
    Hv = h*Hv;
    H = diag(Hv);
    HI = diag(1./Hv);

    % Boundary operators
    e_l = sym(zeros(m,1));
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sym(zeros(m,1));
    d1_l(1:3) = 1/h*[-3/2 2 -1/2];
    d1_r = -rot90(d1_l, 2);

    stencil = [-1/2, 0, 1/2];
    Q = stripeMatrix(stencil, -1:1, m);

    D1 = HI*(Q - 1/2*e_l*e_l' + 1/2*e_r*e_r');

    % Second derivative, 1st order accurate at first boundary points
    M = sym(zeros(m,m));


    c = sym('c', [m,1], 'real');


    scheme_width = 3;
    scheme_radius = (scheme_width-1)/2;
    r = (1+scheme_radius):(m-scheme_radius);


    Mm1 = -c(r-1)/2 - c(r)/2;
    M0  =  c(r-1)/2 + c(r)   + c(r+1)/2;
    Mp1 =            -c(r)/2 - c(r+1)/2;

    M(r,:) = diags([Mm1 M0 Mp1],0:2*scheme_radius,length(r),m);

    M(1:2,1:2) = [c(1)/2 + c(2)/2 -c(1)/2 - c(2)/2; -c(1)/2 - c(2)/2 c(1)/2 + c(2) + c(3)/2;];
    M(m-1:m,m-1:m) = [c(m-2)/2 + c(m-1) + c(m)/2 -c(m-1)/2 - c(m)/2; -c(m-1)/2 - c(m)/2 c(m-1)/2 + c(m)/2;];
    M = 1/h*M;

    D2 = HI*(-M - c(1)*e_l*d1_l' + c(m)*e_r*d1_r');

    r = struct();
    r.Hv = Hv(1);
    r.e_l = e_l(1:3);
    r.d1_l = d1_l(1:3);
    r.D1_inner = stencil;
    r.D1_left_boundary = D1(1,:);
    r.D1_right_boundary = D1(end,:);
    r.D2 = D2;
    r.D2_boundary = D2(1:2,:);
    r.c = c;
    r.h = h;
end
