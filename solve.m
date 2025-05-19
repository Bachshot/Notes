classdef solve
    methods(Static)
        function sol = grow(par)
            sol = struct();

            %% Parameters
            beta = par.beta;
            alpha = par.alpha;
            delta = par.delta;
            sigma = par.sigma;
            gamma = par.gamma;
            nu = par.nu;
            slen = par.slen;
            klen = par.klen;
            Alen = par.Alen;

            kgrid = par.kgrid;
            Agrid = par.Agrid;
            pmat = par.pmat;

            kmat = repmat(kgrid, 1, Alen);
            Amat = repmat(Agrid, klen, 1);

            %% Preallocate
            n0 = zeros(klen, klen, Alen, slen);
            v0 = zeros(klen, Alen, slen); % TEMPORARY INIT for test
            v1 = zeros(klen, Alen, slen);
            k1 = zeros(klen, Alen, slen);
            n1 = zeros(klen, Alen, slen);

            %% Labor supply
            fprintf('------------Solving for Labor Supply.------------\n\n')

            opts = optimset('TolX', 1e-6, 'Display', 'off');  %
            
            for h1 = 1:klen
                for h2 = 1:klen
                    for h3 = 1:Alen
                        for s = 1:par.slen
                            k = kgrid(h1);
                            kp = kgrid(h2);
                            A = Agrid(h3);
                            w = par.w(s);
                            i = kp - (1 - delta) * k;
            
                            fn = @(n) ...
                                -model.utility(w * A * k^alpha * n^(1 - alpha) - ...
                                par.lambda * (w * A * k^alpha * n^(1 - alpha))^(1 - par.tau) - i, ...
                                n, par, 1);  %% g = 1 to decouple it
            
                            n0(h1,h2,h3,s) = fminbnd(fn, 0.001, 0.99, opts);  %% tighten domain
                        end
                    end
                end
            end
            fprintf('------------Labor Supply Done.------------\n\n')

            %% VFI loop
            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;

            fprintf('------------Beginning Value Function Iteration.------------\n\n')

            while diff > crit && iter < maxiter
                for p = 1:klen
                    for j = 1:Alen
                        for s = 1:slen
                            nvec = squeeze(n0(p,:,j,s))';
                            y = par.w(s) * Agrid(j) * (kgrid(p)^alpha) * (nvec.^(1 - alpha));
                            kprime = kgrid;
                            i = kprime - (1 - delta) * kgrid(p);
                            T = par.lambda * y.^(1 - par.tau);
                            c = y - T - i;
                            g = max(mean(T(:)), 1e-4);

                            feasible = (c > 1e-6) & (i >= 0);
                            vall = -Inf(length(kgrid), 1); %% ✅ required init

                            if any(feasible)
                                cvec = c(feasible);
                                nvec_ok = nvec(feasible);
                                vnext = squeeze(v0(:,:,s)) * pmat(j,:)';
                                vall(feasible) = model.utility(cvec(:), nvec_ok(:), par, g) + beta * vnext(feasible);
                            end

                            [vmax, ind] = max(vall);
                            v1(p,j,s) = vmax;
                            k1(p,j,s) = kgrid(ind);
                            n1(p,j,s) = n0(p,ind,j,s);

                            if iter == 0 && p == 1 && j == 1 && s == 1
                                fprintf('w: %.2f | y: %.4f | T: %.4f | i: %.4f | c: %.4f | g: %.4f\n', ...
                                    par.w(s), y(1), T(1), i(1), c(1), g);
                                fprintf('Valid c: %d of %d\n', sum(feasible), length(feasible));
                                fprintf('vmax: %.6f\n', vmax);
                            end
                        end
                    end
                end

                diff = norm(v1(:) - v0(:));
                v0 = v1;
                iter = iter + 1;

                fprintf('Iteration %d: diff = %.10f\n', iter, diff)
            end

            fprintf('\nConverged in %d iterations.\n', iter)
            fprintf('------------End of Value Function Iteration.------------\n')

            %% Store results
            sol.y = zeros(klen, Alen, slen);
            for s = 1:slen
                sol.y(:,:,s) = par.w(s) * Amat .* (kmat.^alpha) .* (n1(:,:,s).^(1 - alpha));
            end
            sol.k = k1;
            sol.n = n1;
            sol.i = k1 - ((1 - delta) * repmat(kgrid, 1, Alen, slen));
            sol.c = sol.y - sol.i;
            sol.v = v1;
        end
    end
end
