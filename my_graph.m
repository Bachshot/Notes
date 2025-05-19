%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim)   
    %% Plot production function for Skill 1 and Skill 2 at A=1


            figure (1)
            hold on
            A_indices = [1, round(par.Alen/2), par.Alen];  % A = low, mid, high
            
            colors = lines(par.slen * length(A_indices));  % Generate enough unique colors
            c = 1;
            for s = 1:par.slen
                for j = A_indices
                    plot(par.kgrid, sol.y(:,j,s), ...
                        'DisplayName', ['Skill ', num2str(s), ', A=', num2str(j)], ...
                        'Color', colors(c,:));
                    c = c + 1;
                end
            end
            
            xlabel({'$k_t$'}, 'Interpreter','latex')
            ylabel({'$y_t$'}, 'Interpreter','latex')
            title('Production Function: Multiple Skills and A States')
            legend show
            hold off

            figure (2)
            hold on
            c = 1;
            for s = 1:par.slen
                for j = A_indices
                    plot(par.kgrid, sol.k(:,j,s), ...
                        'DisplayName', ['Skill ', num2str(s), ', A=', num2str(j)], ...
                        'Color', colors(c,:));
                    c = c + 1;
                end
            end
            xlabel({'$k_t$'}, 'Interpreter','latex')
            ylabel({'$k_{t+1}$'}, 'Interpreter','latex')
            title('Capital Policy Function')
            legend show
            hold off

            figure (3)
            hold on
            c = 1;
            for s = 1:par.slen
                for j = A_indices
                    plot(par.kgrid, sol.c(:,j,s), ...
                        'DisplayName', ['Skill ', num2str(s), ', A=', num2str(j)], ...
                        'Color', colors(c,:));
                    c = c + 1;
                end
            end
            xlabel({'$k_t$'}, 'Interpreter','latex')
            ylabel({'$c_t$'}, 'Interpreter','latex')
            title('Consumption Policy Function')
            legend show
            hold off

            figure (4)
            hold on
            c = 1;
            for s = 1:par.slen
                for j = A_indices
                    plot(par.kgrid, sol.i(:,j,s), ...
                        'DisplayName', ['Skill ', num2str(s), ', A=', num2str(j)], ...
                        'Color', colors(c,:));
                    c = c + 1;
                end
            end
            xlabel({'$k_t$'}, 'Interpreter','latex')
            ylabel({'$i_t$'}, 'Interpreter','latex')
            title('Investment Policy Function')
            legend show
            hold off

            figure(5)
            hold on
            c = 1;
            for s = 1:par.slen
                for j = A_indices
                    plot(par.kgrid, sol.n(:,j,s), ...
                        'DisplayName', ['Skill ', num2str(s), ', A=', num2str(j)], ...
                        'Color', colors(c,:));
                    c = c + 1;
                end
            end
            xlabel({'$k_t$'}, 'Interpreter','latex')
            ylabel({'$n_t$'}, 'Interpreter','latex')
            title('Labor Supply Policy Function')
            legend show
            hold off

            figure (6)
            hold on
            c = 1;
            for s = 1:par.slen
                for j = A_indices
                    plot(par.kgrid, sol.v(:,j,s), ...
                        'DisplayName', ['Skill ', num2str(s), ', A=', num2str(j)], ...
                        'Color', colors(c,:));
                    c = c + 1;
                end
            end
            xlabel({'$k_t$'}, 'Interpreter','latex')
            ylabel({'$v(k_t, A_t)$'}, 'Interpreter','latex')
            title('Value Function')
            legend show
            hold off
        end
    end
end
       