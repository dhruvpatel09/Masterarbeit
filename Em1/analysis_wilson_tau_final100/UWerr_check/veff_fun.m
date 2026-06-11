function y = veff_fun(A, delta_tau)
% A(1) = <W(tau)>
% A(2) = <W(tau_next)>
% Veff = log(W(tau) / W(tau_next)) / delta_tau

y = log(A(1) / A(2)) / delta_tau;
end