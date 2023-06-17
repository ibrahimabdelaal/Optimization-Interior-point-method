function [X_values,Fmin_values,X_S]=Central_path(f,A,b,tolerance,fig)
%%This implementaion solves the whole system in order to get deltas
%%it solves Matrix*Jacpoin=[-rc,-rb,-rxs]
%%Not the best in terms of compution
sa=size(A);
m=sa(1);
A=[A eye(m,m)];  
f=[f zeros(1,m)];
sa=size(A);
Fmin_values=[];
X_values=[];
X_S=[];
%%%Itilaiziatio for  x , s ,y
sigma=.99;
alpha=.5;
[x,s,y]=Intilize_points(f,A,b);
Meo=(x*s')/sa(2);

e=ones(1,sa(2));

for i=1:10000
    Meo=(x*s')/sa(2);
   
    if (Meo)<tolerance ,break ,end

    
    Dx=diag(x);
    Ds=diag(s);
    rxs=(Dx*Ds*e'-(sigma*Meo*e'));
    rc=A'*y'+s'-f';    %dim(3,1)
    rb=A*x'-b';
    siz=size(Dx);
    L=[zeros(sa(2),sa(2)), A' ,eye(sa(2),siz(1));
        A  , zeros(sa(1),sa(1)), zeros(sa(1),siz(1)) ;
        Ds  ,zeros(siz(1),sa(1)),  Dx ];

    h=[-rc;-rb;-rxs];
    delta=L\h;
    deltax=delta(1:sa(2));
    deltay=delta(sa(2)+1:sa(2)+sa(1));
    deltas=delta(end-sa(2)+1:end);
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
    title('Central_path_fixed-central path-feasible region 2d only')
    xlabel('X1')
    ylabel('X2')
    figure()
    p=plot(1:1:length(Fmin_values),Fmin_values);
    p.Marker = '*';
    title('Central_path_fixed-Objective function reduction vs iteration ')
    ylabel('Objective function')
    xlabel('iteration')
    figure()
    p=plot(X_S(:,1),X_S(:,2));
    p.Marker = '*';
    title('Central_path_fixed-complemorty condition ')
    xlabel('X1S1')
    ylabel('X2S2')
end

end