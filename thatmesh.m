f=@(x) 0.75*exp(-1/4*(9*x(1)-2)^2-1/4*(9*x(2)-2)^2)+ ...
    0.75*exp(-(9*x(1)+1)^2/49-(9*x(2)+1)/10)+ ...
    0.5*exp(-(9*x(1)-7)^2/4-(9*x(2)-3)^2/4)- ...
    0.2*exp(-(9*x(1)-4)^2-(9*x(2)-7)^2);

t_evaluate=0:0.01:1;
[X,Y]=meshgrid(t_evaluate);

z=zeros(size(X));

for i=1:length(t_evaluate)
    for j=1:length(t_evaluate)
        Z(i,j)=f([X(i,j);Y(i,j)]);
    end
end
figure

surf(X,Y,Z);
hold on

%% 
% sample sites
HowManySamples=1600;

% H=haltonseq(HowManySamples,2)';
Hq=haltonset(2);
H=net(Hq,HowManySamples)';


Z_samples=zeros(HowManySamples,1);

for i=1:HowManySamples
    Z_samples(i)=f(H(:,i));
end

plot3(H(1,:),H(2,:),Z_samples,'o','LineWidth',6);
% hist(Z_samples,32);
%%
% kernel function
mu=0.05;
K=@(x,y)exp(-1/mu*norm(x-y)^2);
% 
% mu=2;
% K=@(x,y)exp(-1/mu*x'*y);


% % Gram Matrix
% tic
% GramMatrix=zeros(HowManySamples);
% 
% for i=1:HowManySamples
%     for j=1:HowManySamples
%         GramMatrix(i,j)=K(H(:,i),H(:,j));
%     end
% end
% toc


tic
GramMatrix2=exp(-1/mu*(repmat(diag(H'*H)',HowManySamples,1) ...
    -2*(H')*H ...
    + repmat(diag(H'*H)',HowManySamples,1)'));
toc

% weights


Weights=pinv(GramMatrix2)*Z_samples;

% evaluating the approximation

BigEval=[X(:),Y(:)]';
Z_approximation=Weights'*exp(-1/mu*(repmat(diag(H'*H),1,size(BigEval,2))-2*H'*BigEval+repmat(diag(BigEval'*BigEval)',HowManySamples,1)));
Z_approximation02=reshape(Z_approximation,length(t_evaluate),[]);        


% tic
% Z_approximation=zeros(size(X));
% 
% for i=1:length(t_evaluate)
%     for j=1:length(t_evaluate)
% %         Z_approximation(i,j)=Weights'*exp(-1/mu*(diag(H'*H)-2*H'*[X(i,j);Y(i,j)]+[X(i,j);Y(i,j)]'*[X(i,j);Y(i,j)]));
%         
%         sum=0;
%         for ii=1:HowManySamples
%             sum=sum+Weights(ii)*K([X(i,j);Y(i,j)],H(:,ii));
%         end
%         Z_approximation(i,j)=sum;
%     end
% end
% toc
% plot

figure
surf(X,Y,Z_approximation02);
