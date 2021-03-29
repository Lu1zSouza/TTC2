#Método da Secante

function [p] = secante(f,a,b,TOL,N);
x0 = a;
x1 = b;
x2 = 0;
erel = TOL +1;
i = 0;
err = [];
printf("n \t  x0 \t\t  f(x) \t\t  x1 \t\t  f(x1) \t x2 \t\terr(%%)\n");
while erel > TOL && i<N
  x2 = (x0*f(x1) - x1*f(x0))/(f(x1) - f(x0));
  if x2!=0
      erel = abs((x2-x1)/x2);
      err(i+1) = erel;
  endif
  i = i +1;
  printf("%d \t %f \t %f \t %f \t %f \t %f\t %f\n", i, x0, f(x0), x1, f(x1), x2, erel*100);
  x0 = x1;
  x1 = x2;
endwhile  
  printf("A resposta é %f\n", x2);
  i=1:i;
    plot(i,err(i));
    xlabel ("Quantidade de Iterações");
    ylabel ("Erro Relativo");
    title ("Evolução do Erro Relativo");

end
