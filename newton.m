%Método de Newton-Raphson

function [p] = newton(f, dfdx, x0, TOL, N)
    x = x0;
    x1 = x0;
    fx = f(x);
    i = 1;
    erel = 0;
    err = [];
    printf("n \t  x0 \t\t  f(x) \t\t  f'(x) \t  x1 \t\terr(%%)\n");
    while abs(fx) > TOL && i <= N
        try
            x = x - (fx)/dfdx(x);
            if x!= 0
              erel = abs((x - x1)/x);
              err(i) = erel;
            endif
            printf("%d \t % f \t % f \t % f \t % f \t%f\n", i, x1, fx, dfdx(x), x, erel*100);
            x1 = x;
        catch
            printf('Erro! - Derivada 0 para x = \n', x)
            exit(1)
        end
        fx = f(x);
        i = i + 1;
    end
    i=1:i-1;
    plot(i,err(i));
    xlabel ("Quantidade de Iterações");
    ylabel ("Erro Relativo");
    title ("Evolução do Erro Relativo");
    if abs(fx) > TOL
        i = -1;
    end
    p = x;
    printf("A resposta é %f\n", p);
    
end
