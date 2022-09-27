        new_num=500;
%      cut_index=randsrc(1,1,randperm(k));
%      cut_point_norm=failed_sample_norm(2,:);
%      cut_point_norm=cut_mean;
%         new_sample_norm_tmp=normrnd(repmat(failed_sample_norm(1,i),new_num,1),repmat(ones(1,1),new_num,1));
        new_sample_norm_tmp=normrnd(repmat(cut_point_norm(1,i),new_num,1),repmat(ones(1,1),new_num,1));
        new_sample_norm=repmat(cut_point_norm,new_num,1);
        new_sample_norm(:,i)=new_sample_norm_tmp;
        new_sample_run=repmat(mu,new_num,1)+new_sample_norm.*repmat(sigma,new_num,1);
        [new_delay]=run_file(analysis.hspicepath,new_sample_run,new_num)-ft;
        tspld=ppval(dmodel_spl,new_sample_norm(:,i));
%         tpoly=polyval(dmodel_pol,new_sample_norm(:,i));
        
        % added 
        tpoly2=polyval(dmodel_pol_2,new_sample_norm(:,i));
        tpoly3=polyval(dmodel_pol_3,new_sample_norm(:,i));
        tpoly4=polyval(dmodel_pol_4,new_sample_norm(:,i));
        tpoly5=polyval(dmodel_pol_5,new_sample_norm(:,i));
        
        tgp=predict(gprMdl,new_sample_norm(:,i));
       %% visualize
        vect=[new_sample_norm(:,i),new_delay+ft];
        vect_sort=sort(vect,1);
%         plot(vect_sort(:,1),vect_sort(:,2),'k*','Marksize',40)
        plot(vect_sort(:,1),vect_sort(:,2),'.','MarkerSize',20,'Color','k');
        hold on
        vecSPL=[new_sample_norm(:,i),tspld+ft];
        vecSPL_sort=sort(vecSPL,1);
        plot(vecSPL_sort(:,1),vecSPL_sort(:,2),'r','linewidth',2)
        hold on
        % added 
        vecPoly2=[new_sample_norm(:,i),tpoly2+ft];
        vecPoly_sort2=sort(vecPoly2,1);
        plot(vecPoly_sort2(:,1),vecPoly_sort2(:,2),'g','linewidth',2)
        hold on
        vecPoly3=[new_sample_norm(:,i),tpoly3+ft];
        vecPoly_sort3=sort(vecPoly3,1);
        plot(vecPoly_sort3(:,1),vecPoly_sort3(:,2),'b','linewidth',2)
        hold on
        vecPoly4=[new_sample_norm(:,i),tpoly4+ft];
        vecPoly_sort4=sort(vecPoly4,1);
        plot(vecPoly_sort4(:,1),vecPoly_sort4(:,2),'m','linewidth',2)
        hold on
        vecPoly5=[new_sample_norm(:,i),tpoly5+ft];
        vecPoly_sort5=sort(vecPoly5,1);
        plot(vecPoly_sort5(:,1),vecPoly_sort5(:,2),'c','linewidth',2)
%         hold on
%         vecGP=[new_sample_norm(:,i),tgp];
%         vecGP_sort=sort(vecGP,1);
%         plot(vecGP_sort(:,1),vecGP_sort(:,2),'y','linewidth',2)
        xlabel('"Vth0" variation (normalized) ')
        ylabel('Read Delay')
%         ylabel('Read Delay','FontWeight','bold')
%         ylabel({$$'f_0'},'Interpreter','latex','FontSize',18)
%             ylabel('\alpha^{2}+\beta_{0}^{2}=1')
% title('$\frac{f(x)}{x}$','interpreter','latex', 'FontSize', 18);
% ylabel('$\f(x)$','interpreter','latex', 'FontSize', 18);

% ylabel({'$f_(x_i,\bf{x_{0}^{i}})-f_{0}$'} ,'interpreter','latex','FontWeight','bold')
%         ylabel('frac{f(x)-f_{0}}','FontWeight','bold');
        set(gca,'FontSize',30)
        grid on
        legend({'HSPICE Results','RBF interpolation','2-order Polynomial','3-order Polynomial','4-order Polynomial','5-order Polynomial'},'FontSize',22);
