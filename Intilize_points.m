function [x,s,y]=Intilize_points(f,A,b)
x_hat=A'*(A*A')^-1*b';
y_hat=(A*A')^-1*A*f';
s_hat=f'-A'*y_hat;

aplha_x=max(-1.5*min(x_hat),0);
aplha_s=max(-1.5*min(s_hat),0);
 

x=x_hat+aplha_x;
s=s_hat+aplha_s;
e=ones(1,length(s));

aplha_hat_x=x'*s/(2*e*s);
aplha_hat_s=x'*s/(2*e*x);

x=x+aplha_hat_x;
s=s+aplha_hat_s;
x=x';
s=s';
y=y_hat';





end