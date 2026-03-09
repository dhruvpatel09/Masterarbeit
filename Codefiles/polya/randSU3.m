function u = randSU3()
   persistent Ta
   if isempty(Ta)
      Ta{1} = [0 1 0; 1 0 0; 0 0 0];
      Ta{2} = [0 -1i 0; 1i 0 0; 0 0 0];
      Ta{3} = [1 0 0; 0 -1 0; 0 0 0];
      Ta{4} = [0 0 1; 0 0 0; 1 0 0];
      Ta{5} = [0 0 -1i; 0 0 0; 1i 0 0];
      Ta{6} = [0 0 0; 0 0 1; 0 1 0];
      Ta{7} = [0 0 0; 0 0 -1i; 0 1i 0];
      Ta{8} = [1 0 0; 0 1 0; 0 0 -2]/sqrt(3);
   end
   r = randn(8,1);
   u = expm(1i*(r(1)*Ta{1}+r(2)*Ta{2}+r(3)*Ta{3}+r(4)*Ta{4}...
               +r(5)*Ta{5}+r(6)*Ta{6}+r(7)*Ta{7}+r(8)*Ta{8}));
end