function [X_values,Fmin_values,X_S]=Mehrotra_Predictor(f,A,b,n,tolerance,fig)
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
if (n<.9 || n>1) ,disp("ERROR : Invalid n shoud be between .9,1"); return ,end
%%%Itilaiziatio for  x , s ,y
e=ones(1,sa(2));
[x,s,y]=Intilize_points(f,A,b);
for i=1:10000 
     Meo=(x*s')/sa(2);
   % if (x*s'/sa(2))<0 , Meo=.01;,end
    if (Meo)<tolerance ,break ,end

  %%%  Prediction : Solvig the affine problem 
    
    Dx=diag(x);
    Ds=diag(s);
    rxs_affine=(Dx*Ds*e');
    rc=A'*y'+s'-f';    %dim(3,1)
    rb=A*x'-b';
    siz=size(Dx);
    L=[zeros(sa(2),sa(2)), A' ,eye(sa(2),siz(1));
        A  , zeros(sa(1),sa(1)), zeros(sa(1),siz(1)) ;
        Ds  ,zeros(siz(1),sa(1)),  Dx ];

    h=[-rc;-rb;-rxs_affine];
    %%% Solvinh L*delta=h
    delta_affine=L\h;
    deltax_affine=delta_affine(1:sa(2));
    deltas_affine=delta_affine(sa(2)+sa(1)+1:end);
    alpha_prime=min(1,min(-1*x(deltax_affine<0))./deltax_affine(deltax_affine<0));
    alpha_dual=min(1,min(-1*s(deltas_affine<0)./(deltas_affine(deltas_affine<0))));
    if(isempty(alpha_prime))  alpha_prime=1; , else  alpha_prime=alpha_prime(1);, end
    if(isempty(alpha_dual))  alpha_dual=1; , else  alpha_dual=alpha_dual(1);, end
    
    Meo_aff=((x'+alpha_prime*deltax_affine)'*(s'+alpha_dual*deltas_affine))/sa(2);
    %%%% Corrector step 
    sigma=(Meo_aff/Meo)^3;
    rxs=(Dx*Ds*e'+(deltax_affine.*deltas_affine)-(sigma*Meo));
    rc=A'*y'+s'-f';    %dim(3,1)
    rb=A*x'-b';
    L=[zeros(sa(2),sa(2)), A' ,eye(sa(2),siz(1));A  , zeros(sa(1),sa(1)), zeros(sa(1),siz(1)) ;Ds  ,zeros(siz(1),sa(1)),  Dx ];
    h=[-rc;-rb;-rxs];
    delta=L\h;
    deltax=delta(1:sa(2));
    deltay=delta(sa(2)+1:sa(2)+sa(1));
    deltas=delta(sa(2)+sa(1)+1:end);   
    alpha_prime=min(-1*x(deltax<0))./deltax(deltax<0);
    alpha_dual=min(-1*s(deltas<0)./(deltas(deltas<0)));
    if(isempty(alpha_prime))  alpha_prime=.9; , else  alpha_prime=min(.9,n*alpha_prime(1));, end
    if(isempty(alpha_dual))  alpha_dual=.9; , else  alpha_dual=min(.9,n*alpha_dual(1));, end
    %%Update step 
    x=x+(alpha_prime*deltax');
    s=s+(alpha_dual*deltas');
    y=y+alpha_dual*deltay';
    
    Fmin_values(i)=f*x';
    X_values(i,:)=x(1:end-m);
    X_S(i,:)=((x.*s));

end
if fig 
    figure()
    A_orig=A(:,1:end-m);
    boundries=b'./A_orig;
    hold on
    %plotting boundries for two variables only
    for j=1:size(A_orig,1)
        u=plot([boundries(j,1),0],[0,boundries(j,2)]);
       
    end
    
    %a=plot(Fmin_values,1:1:length(Fmin_values));
    p=plot(X_values(:,1),X_values(:,2));
    p.Marker = '*';
    title('Mehrotra_Predictor-central path-contraints 2d only')
    xlabel('X1')
    ylabel('X2')
    legend(p,'central path') 
    legend(u,'contraints') 
    figure()
    p=plot(1:1:length(Fmin_values),Fmin_values);
    p.Marker = '*';
    title('Mehrotra_Predictor-Objective function reduction vs iteration ')
    ylabel('Objective function')
    xlabel('iteration')
    figure()
    
    p=plot(X_S(:,1),X_S(:,2));
    p.Marker = '*';
    title('Mehrotra_Predictor-complemorty condition ')
    xlabel('X1S1')
    ylabel('X2S2')
end

end