function a=amn(m,n);
    if m==n
        a = -n^2*pi^2/2;
    else
        a = - n^2 * ( ((-1)^(m-n) -1)/(m-n)^2 - ((-1)^(m+n) -1)/(m+n)^2 );
    end
