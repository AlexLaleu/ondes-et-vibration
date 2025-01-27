function b=bmn(m,n);
    if m==n
        b = 0 ;
    else
        b = - n * ( ((-1)^(m-n) -1)/(m-n) + ((-1)^(m+n) -1)/(m+n) );
    end
