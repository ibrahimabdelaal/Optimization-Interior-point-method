function [X_values,Fmin_values,X_S]=Mehrotra_Predictor_2(f,A,b,n,tolerance,fig)
%%This implementaion solves the only the normal equation to get deltay
%%then calculate deltas from delta y 
%%then calculate deltax from delta s
%%better than Mehrotra_Predictor in terms of effciency
sa=size(A);
num_var=sa(2);
m=sa(1);
A=[A eye(m,m)];  
f=[f zeros(1,m)];
sa=size(A);
Fmin_values=[];
X_values=[];
X_S=[];

%%%Intilaiziatio for  x , s ,y
[x,s,y]=Intilize_points(f,A,b);
e=ones(1,sa(2));

for i=1:5000
    if (x*s'/sa(2))<tolerance ,break ,end
    Dx=diag(x);
    Ds=diag(s);
    rxs_affine=(Dx*Ds*e');
    rc=A'*y'+s'-f';    %dim(3,1)
    rb=A*x'-b';
    inv_s=inv(Ds);
    D2=inv_s*Dx;
    deltay_affine=(A*D2*A')\((-rb-A*Ds^-1*Dx*rc)+(A*Ds^-1*rxs_affine));
    deltas_affine=-rc-A'*deltay_affine;
    deltax_affine=(-Ds^-1*rxs_affine)-(Ds^-1*Dx*deltas_affine);
    alpha_prime=min(1,min(-1*x(deltax_affine<0))./deltax_affine(deltax_affine<0));
    alpha_dual=min(1,min(-1*s(deltas_affine<0)./(deltas_affine(deltas_affine<0))));
    if(isempty(alpha_prime))  alpha_prime=1; , else  alpha_prime=alpha_prime(1);, end
    if(isempty(alpha_dual))  alpha_dual=1; , else  alpha_dual=alpha_dual(1);, end
    Meo_aff=((x+alpha_prime*deltax_affine')*(s+alpha_dual*deltas_affine')')/sa(2);
    Meo=(x*s')/sa(2);
    sigma=(Meo_aff/Meo)^3;
    rxs=(Dx*Ds*e'+(deltax_affine.*deltas_affine)-(sigma*Meo*e'));
    rc=A'*y'+s'-f';    %dim(3,1)
    rb=A*x'-b';
    deltay=(A*D2*A')\((-rb-A*Ds^-1*Dx*rc)+(A*Ds^-1*rxs));
    deltas=-rc-A'*deltay;
    deltax=(-inv_s*rxs)-(inv_s*Dx*deltas);
    alpha_prime=(min(-1*x(deltax<0))./deltax(deltax<0));
    alpha_dual=(min(-1*s(deltas<0)./(deltas(deltas<0))));
    if(isempty(alpha_prime))  alpha_prime=1; , else  alpha_prime=min(.9,n*alpha_prime(1));, end
    if(isempty(alpha_dual))  alpha_dual=1; , else  alpha_dual=min(.9,n*alpha_dual(1));, end
    
   
    x=x+(alpha_prime*deltax');
    s=s+(alpha_dual*deltas');
    y=y+alpha_dual*deltay';
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
        u=plot([boundries(j,1),0],[0,boundries(j,2)]);
       
    end
    
    %a=plot(Fmin_values,1:1:length(Fmin_values));
    p=plot(X_values(:,1),X_values(:,2));
    p.Marker = '*';
    title('Mehrotra_Predictor2-central path-contraints 2d only')
    xlabel('X1')
    ylabel('X2')
    legend(p,'central path') 
    legend(u,'contraints') 
    figure()
    p=plot(1:1:length(Fmin_values),Fmin_values);
    p.Marker = '*';
    title('Mehrotra_Predictor2-Objective function reduction vs iteration ')
    ylabel('Objective function')
    xlabel('iteration')
    figure()
    
    p=plot(X_S(:,1),X_S(:,2));
    p.Marker = '*';
    title('Mehrotra_Predictor2-complemorty condition ')
    xlabel('X1S1')
    ylabel('X2S2')

end

end