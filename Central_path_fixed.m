function [X_values,Fmin_values,X_S]=Central_path_fixed(f,A,b,tolerance,fig)
%%%This implementaion solves the only the normal equation to get deltay
%%then calculate deltas from delta y 
%%then calculate deltax from delta s
%%better than central_path function in terms of effciency
sa=size(A);
m=sa(1);
A=[A eye(m,m)];  
f=[f zeros(1,m)];
sa=size(A);
Fmin_values=[];
X_values=[];
X_S=[];
%%%Itilaiziatio for  x , s ,y
sigma=.8;
alpha=.01;
[x,s]=Intilize_points(f,A,b);

Meo=(x*s')/sa(2);
y=ones(1,sa(1));
e=ones(1,sa(2));

for i=1:10000
    if (x*s'/sa(2))<tolerance ,break ,end
    Meo=(x*s')/sa(2);
    Dx=diag(x);
    Ds=diag(s);
    rxs=(Dx*Ds*e'-(sigma*Meo*e'));
    rc=A'*y'+s'-f';    %dim(3,1)
    rb=A*x'-b';
    inv_s=inv(Ds);
    D2=inv_s*Dx;
    deltay=A*D2*A'\((-rb-A*inv_s*Dx*rc)+(A*inv_s*rxs));
    %%%delta s
    deltas=-rc-A'*deltay;
    deltax=(-inv_s*rxs)-(inv_s*Dx*deltas);
    %%% updatig 
    x=x+(alpha*deltax');
    s=s+(alpha*deltas');
    y=y+alpha*deltay';
    Fmin_values(i)=f*x';
    X_values(i,:)=x(1:end-m);
    X_S(i,:)=(x.*s);
 
end
if fig
    figure()
    A_orig=A(:,1:end-m);
    boundries=b'./A_orig;
    hold on
    %plotting boundries for two variables only
    for j=1:size(A_orig,1)
        plot([boundries(j,1),0],[0,boundries(j,2)]);
       
    end
    
    %a=plot(Fmin_values,1:1:length(Fmin_values));
    p=plot(X_values(:,1),X_values(:,2));
    p.Marker = '*';
    title('Central_path_fixed_lesscompution-central path-feasible region 2d only')
    xlabel('X1')
    ylabel('X2')
    figure()
    p=plot(1:1:length(Fmin_values),Fmin_values);
    p.Marker = '*';
    title('Central_path_fixed_lesscompution-Objective function reduction vs iteration ')
    ylabel('Objective function')
    xlabel('iteration')
    figure()
    p=plot(X_S(:,1),X_S(:,2));
    p.Marker = '*';
    title('Central_path_fixed_lesscompution-complemorty condition ')
    xlabel('X1S1')
    ylabel('X2S2')




end



end