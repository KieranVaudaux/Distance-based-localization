%% Dimension 3
clc;
clear;
num_exp = 12;
n = 35;
tic;
N = [n];
d = 6;
M = [10,12,16,18,20,24,26];
Ratio = [0.25,0.28,0.33,0.37,0.43,0.49,0.57];

r_minus_d_list = [0];

for i=1:length(N)
    n = N(i);
    for j=1:length(r_minus_d_list)
        

        r_minus_d = r_minus_d_list(j);

        % Graph 1
        cost_relax_list1 = zeros(num_exp,length(M(i,:)));
        cost_proj_list1 = zeros(num_exp,length(M(i,:)));
        m_list1 = zeros(num_exp,length(M(i,:)));

        for k=1:length(M(i,:))
           

            for l=1:num_exp
                 m = M(i,k);
                [X, G, C, D, edges, A] = graphSphere(n,m,d,true,true);
                m = size(D,1)
                m_list1(l,k) = m;
                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
                
                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 60;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list1(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list1(l,k) = abs(cost(Y_proj));
                end

            end
        end

        % Graph 2
        cost_relax_list2 = zeros(num_exp,length(M(i,:)));
        cost_proj_list2 = zeros(num_exp,length(M(i,:)));
        m_list2 = zeros(num_exp,length(Ratio));

        for k=1:length(M(i,:))
            m = M(i,k);

            for l=1:num_exp
                 m = M(i,k);
                [X, G, C, D, edges, A] = graphSphere(n,m,d,false,true);
                m = size(D,1);
                m_list2(l,k) = m;
                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');

                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 60;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list2(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list2(l,k) = abs(cost(Y_proj));
                end

                % Infos for convergence of singular value

                
            end
        end

        % Graph 3
        cost_relax_list3 = zeros(num_exp,length(M(i,:)));
        cost_proj_list3 = zeros(num_exp,length(M(i,:)));
        m_list3 = zeros(num_exp,length(Ratio));

        for k=1:length(M(i,:))
            m = M(i,k);

            for l=1:num_exp
                 m = M(i,k);
                [X, G, C, D, edges, A] = graphSphere(n,m,d,true,false);
                m = size(D,1);
                m_list3(l,k) = m;
                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');

                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 60;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list3(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list3(l,k) = abs(cost(Y_proj));
                end

                
            end
        end

         % Graph 4
        cost_relax_list4 = zeros(num_exp,length(M(i,:)));
        cost_proj_list4 = zeros(num_exp,length(M(i,:)));
        m_list4 = zeros(num_exp,length(Ratio));

        for k=1:length(M(i,:))
            ratio = Ratio(i,k);

            for l=1:num_exp
                X = normrnd(0,1,n,d);
                [X, G, C, D, edges, A] = randomGraph(X,ratio);
                m = size(D,1);
                m_list4(l,k) = m;

                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');

                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 30;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list4(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list4(l,k) = abs(cost(Y_proj));
                end

            end
        end


          
        if r_minus_d > 0
            clear fig;
            [mean_,ind] = sort(mean(m_list1,1));
            mean_cost = mean(cost_proj_list1,1);
 
            fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the ball');
            hold on;
            
           [mean_,ind] = sort(mean(m_list2,1));
           mean_cost = mean(cost_proj_list2,1);
           semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the sphere');
           
           [mean_,ind] = sort(mean(m_list3,1));
            mean_cost = mean(cost_proj_list3,1);
            semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the ball');
            
            [mean_,ind] = sort(mean(m_list4,1));
            mean_cost = mean(cost_proj_list4,1);
            semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');
            
            legend();
            xlabel('Number of edges');
            ylabel('Cost');
            hold off;
            name = strcat('FigData/DimensionD/mean_cost_proj_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png')
            saveas(fig,name);
            
        end
    
        
        clear fig;
        [mean_,ind] = sort(mean(m_list1,1));
        mean_cost = mean(cost_relax_list1,1);
        fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the ball');
        hold on;
        
       [mean_,ind] = sort(mean(m_list2,1));
       mean_cost = mean(cost_relax_list2,1);
       semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the sphere');
       
       [mean_,ind] = sort(mean(m_list3,1));
        mean_cost = mean(cost_relax_list3,1);
        semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the ball');
        
        [mean_,ind] = sort(mean(m_list4,1));
        mean_cost = mean(cost_relax_list4,1);
        semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');

        legend();
        xlabel('Number of edges');
        ylabel('Cost');
        name = strcat('FigData/DimensionD/mean_cost_relax_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png')
        hold off;
        saveas(fig,name);
        

        if r_minus_d > 0
            clear fig;
            
            [mean_,ind] = sort(mean(m_list1,1));
            mean_cost = median(cost_proj_list1,1);
            fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the ball');
            hold on;
            
           [mean_,ind] = sort(mean(m_list2,1));
           mean_cost = median(cost_proj_list2,1);
           semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the spehre');
           
           [mean_,ind] = sort(mean(m_list3,1));
            mean_cost = median(cost_proj_list3,1);
            semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the ball');
            
            [mean_,ind] = sort(mean(m_list4,1));
            mean_cost = median(cost_proj_list4,1);
            semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');


            legend();
            xlabel('Number of edges');
            ylabel('Cost');
            hold off;
            name = strcat('FigData/DimensionD/median_cost_proj_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png');
            saveas(fig,name);
            
        end
    
        
        clear fig;
        [mean_,ind] = sort(mean(m_list1,1));
        mean_cost = median(cost_relax_list1,1);
        fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the ball');
        hold on;
        
       [mean_,ind] = sort(mean(m_list2,1));
       mean_cost = median(cost_relax_list2,1);
       semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the sphere');
       
       [mean_,ind] = sort(mean(m_list3,1));
        mean_cost = median(cost_relax_list3,1);
        semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the ball');
        
        [mean_,ind] = sort(mean(m_list4,1));
        mean_cost = median(cost_relax_list4,1);
        semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');

        
        legend();
        xlabel('Number of edges');
        ylabel('Cost');
        hold off;
        name = strcat('FigData/DimensionD/median_cost_relax_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png');
        saveas(fig,name);
        

    end
end
toc


%% Dimension 3
clc;
clear;
num_exp = 12;
n = 35;
tic;
N = [n];
d = 6;
M = [10,12,16,18,20,24,26];
Ratio = [0.25,0.28,0.33,0.37,0.43,0.49,0.57];

r_minus_d_list = [3];

for i=1:length(N)
    n = N(i);
    for j=1:length(r_minus_d_list)
        

        r_minus_d = r_minus_d_list(j);

        % Graph 1
        cost_relax_list1 = zeros(num_exp,length(M(i,:)));
        cost_proj_list1 = zeros(num_exp,length(M(i,:)));
        m_list1 = zeros(num_exp,length(M(i,:)));

        for k=1:length(M(i,:))
           

            for l=1:num_exp
                 m = M(i,k);
                [X, G, C, D, edges, A] = graphSphere(n,m,d,true,true);
                m = size(D,1)
                m_list1(l,k) = m;
                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
                
                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 60;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list1(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list1(l,k) = abs(cost(Y_proj));
                end

            end
        end

        % Graph 2
        cost_relax_list2 = zeros(num_exp,length(M(i,:)));
        cost_proj_list2 = zeros(num_exp,length(M(i,:)));
        m_list2 = zeros(num_exp,length(Ratio));

        for k=1:length(M(i,:))
            m = M(i,k);

            for l=1:num_exp
                 m = M(i,k);
                [X, G, C, D, edges, A] = graphSphere(n,m,d,false,true);
                m = size(D,1);
                m_list2(l,k) = m;
                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');

                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 60;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list2(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list2(l,k) = abs(cost(Y_proj));
                end

                % Infos for convergence of singular value

                
            end
        end

        % Graph 3
        cost_relax_list3 = zeros(num_exp,length(M(i,:)));
        cost_proj_list3 = zeros(num_exp,length(M(i,:)));
        m_list3 = zeros(num_exp,length(Ratio));

        for k=1:length(M(i,:))
            m = M(i,k);

            for l=1:num_exp
                 m = M(i,k);
                [X, G, C, D, edges, A] = graphSphere(n,m,d,true,false);
                m = size(D,1);
                m_list3(l,k) = m;
                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');

                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 60;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list3(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list3(l,k) = abs(cost(Y_proj));
                end

                
            end
        end

         % Graph 4
        cost_relax_list4 = zeros(num_exp,length(M(i,:)));
        cost_proj_list4 = zeros(num_exp,length(M(i,:)));
        m_list4 = zeros(num_exp,length(Ratio));

        for k=1:length(M(i,:))
            ratio = Ratio(i,k);

            for l=1:num_exp
                X = normrnd(0,1,n,d);
                [X, G, C, D, edges, A] = randomGraph(X,ratio);
                m = size(D,1);
                m_list4(l,k) = m;

                W = eye(m);

                Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
                Q1 = D*W*D;
                cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');

                manifold = obliquefactory(d+r_minus_d, m,true);
                Y_0 = manifold.rand();
                init = true;
                option.tolgradnorm = 10^-6;
                option.maxtime = 30;
                option.verbosity = false;
                
                [x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);
                cost_relax_list4(l,k) = abs(xcost);

                if r_minus_d > 0
                    [U,S,V] = svd(Y);
                    Y_proj = normr(U*S(:,1:d));
                    cost_proj_list4(l,k) = abs(cost(Y_proj));
                end

            end
        end


         
        if r_minus_d > 0
            clear fig;
            [mean_,ind] = sort(mean(m_list1,1));
            mean_cost = mean(cost_proj_list1,1);
 
            fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the disk');
            hold on;
            
           [mean_,ind] = sort(mean(m_list2,1));
           mean_cost = mean(cost_proj_list2,1);
           semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the circle');
           
           [mean_,ind] = sort(mean(m_list3,1));
            mean_cost = mean(cost_proj_list3,1);
            semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the disk');
            
            [mean_,ind] = sort(mean(m_list4,1));
            mean_cost = mean(cost_proj_list4,1);
            semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');
            
            legend();
            xlabel('Number of edges');
            ylabel('Cost');
            hold off;
            name = strcat('FigData/DimensionD/mean_cost_proj_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png')
            saveas(fig,name);
            
        end
    
        
        clear fig;
        [mean_,ind] = sort(mean(m_list1,1));
        mean_cost = mean(cost_relax_list1,1);
        fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the disk');
        hold on;
        
       [mean_,ind] = sort(mean(m_list2,1));
       mean_cost = mean(cost_relax_list2,1);
       semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the circle');
       
       [mean_,ind] = sort(mean(m_list3,1));
        mean_cost = mean(cost_relax_list3,1);
        semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the disk');
        
        [mean_,ind] = sort(mean(m_list4,1));
        mean_cost = mean(cost_relax_list4,1);
        semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');

        legend();
        xlabel('Number of edges');
        ylabel('Cost');
        name = strcat('FigData/DimensionD/mean_cost_relax_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png')
        hold off;
        saveas(fig,name);
        

        if r_minus_d > 0
            clear fig;
            
            [mean_,ind] = sort(mean(m_list1,1));
            mean_cost = median(cost_proj_list1,1);
            fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the disk');
            hold on;
            
           [mean_,ind] = sort(mean(m_list2,1));
           mean_cost = median(cost_proj_list2,1);
           semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the circle');
           
           [mean_,ind] = sort(mean(m_list3,1));
            mean_cost = median(cost_proj_list3,1);
            semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the disk');
            
            [mean_,ind] = sort(mean(m_list4,1));
            mean_cost = median(cost_proj_list4,1);
            semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');


            legend();
            xlabel('Number of edges');
            ylabel('Cost');
            hold off;
            name = strcat('FigData/DimensionD/median_cost_proj_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png');
            saveas(fig,name);
            
        end
    
        
        clear fig;
        [mean_,ind] = sort(mean(m_list1,1));
        mean_cost = median(cost_relax_list1,1);
        fig = semilogy(mean_,mean_cost(ind),'-square','DisplayName','Local structure on the disk');
        hold on;
        
       [mean_,ind] = sort(mean(m_list2,1));
       mean_cost = median(cost_relax_list2,1);
       semilogy(mean_,mean_cost(ind),'-o','DisplayName','Local structure on the circle');
       
       [mean_,ind] = sort(mean(m_list3,1));
        mean_cost = median(cost_relax_list3,1);
        semilogy(mean_,mean_cost(ind),'-hexagram','DisplayName','Non local structure on the disk');
        
        [mean_,ind] = sort(mean(m_list4,1));
        mean_cost = median(cost_relax_list4,1);
        semilogy(mean_,mean_cost(ind),'-diamond','DisplayName','Random graph');

        
        legend();
        xlabel('Number of edges');
        ylabel('Cost');
        hold off;
        name = strcat('FigData/DimensionD/median_cost_relax_n_',int2str(n),'_rmd_',int2str(r_minus_d),'.png');
        saveas(fig,name);
        

    end
end
toc
