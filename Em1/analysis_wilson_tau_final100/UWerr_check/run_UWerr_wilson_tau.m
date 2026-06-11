clear;
close all;

addpath(pwd);

diary('UWerr_results.txt');
diary on;

fprintf('\n============================================\n');
fprintf('UWerr cross-check for Wilson-tau analysis\n');
fprintf('============================================\n');

raw = dlmread('wilson_tau_per_config.csv', ',', 1, 0);

% CSV columns:
% 1 tau
% 2 nc
% 3 corr_id
% 4 ReW
% 5 ImW
% 6 AbsW

primary_taus = [5, 8];
derived_pairs = [5, 6; 8, 9];

% We run both:
% Stau = 1.5  : UWerr default-like choice
% Stau = 2.0  : closer to our pyerrors script setting
Stau_list = [1.5, 2.0];

for is = 1:length(Stau_list)

    Stau = Stau_list(is);

    fprintf('\n\n##################################################\n');
    fprintf('Running UWerr with Stau = %.1f\n', Stau);
    fprintf('##################################################\n');

    % ------------------------------------------------------------
    % Primary observables ReW(tau)
    % ------------------------------------------------------------
    for k = 1:length(primary_taus)

        tau = primary_taus(k);

        rows = raw(raw(:,1) == tau, :);
        rows = sortrows(rows, 2);

        D = rows(:,4);

        fprintf('\n============================================\n');
        fprintf('Primary observable: ReW(tau=%d)\n', tau);
        fprintf('Ncfg = %d\n', length(D));
        fprintf('============================================\n');

        % Name=0 suppresses UWerr plots, useful for no-X11 Octave.
        [value,dvalue,ddvalue,tauint,dtauint,Qval,Wopt] = ...
            UWerr(D, Stau, [], 0, 1);

        fprintf('value   = %.15e\n', value);
        fprintf('dvalue  = %.15e\n', dvalue);
        fprintf('ddvalue = %.15e\n', ddvalue);
        fprintf('tauint  = %.15e\n', tauint);
        fprintf('dtauint = %.15e\n', dtauint);
        fprintf('Wopt    = %d\n', Wopt);
    end


    % ------------------------------------------------------------
    % Derived observables Veff(tau -> tau_next)
    % ------------------------------------------------------------
    for k = 1:size(derived_pairs,1)

        tau1 = derived_pairs(k,1);
        tau2 = derived_pairs(k,2);
        delta_tau = tau2 - tau1;

        rows1 = raw(raw(:,1) == tau1, :);
        rows2 = raw(raw(:,1) == tau2, :);

        rows1 = sortrows(rows1, 2);
        rows2 = sortrows(rows2, 2);

        if length(rows1(:,2)) ~= length(rows2(:,2))
            error('Different number of configurations for tau pair');
        end

        if any(rows1(:,2) ~= rows2(:,2))
            error('Configuration numbers do not match for tau pair');
        end

        DataPair = [rows1(:,4), rows2(:,4)];

        fprintf('\n============================================\n');
        fprintf('Derived observable: Veff(%d -> %d)\n', tau1, tau2);
        fprintf('Ncfg = %d\n', size(DataPair,1));
        fprintf('============================================\n');

        [value,dvalue,ddvalue,tauint,dtauint,Qval,Wopt] = ...
            UWerr(DataPair, Stau, [], 0, @veff_fun, delta_tau);

        fprintf('value   = %.15e\n', value);
        fprintf('dvalue  = %.15e\n', dvalue);
        fprintf('ddvalue = %.15e\n', ddvalue);
        fprintf('tauint  = %.15e\n', tauint);
        fprintf('dtauint = %.15e\n', dtauint);
        fprintf('Wopt    = %d\n', Wopt);
    end
end

diary off;