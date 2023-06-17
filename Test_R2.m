%%Test cetrntral path
%%General instruction 
%%%please note that input must be in the form 
%%% Min f (X) = c1x1 + c2x2 + · · · + cnxn     (1)
%%% Subject Ax <= b   //Do not add slack variables 
%%%please,Do not add slack variables :the algorithm takes care of it
clear
clc
close all
f=[-1.1,-1];
A=[1,1];
b=6;
n=.91;
sigma=.8;
tolerance=.00001;
fig=0;             %%flag to display figures : set 1 to display
[X_values,Fmin_values,X_S]=Mehrotra_Predictor(f,A,b,n,tolerance,fig);
disp("########### Moherete predictor results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Mehrotra_Predictor_2(f,A,b,n,tolerance,fig);
disp("########### Moherete 2 (better compution effciency predictor results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))


[X_values,Fmin_values,X_S]=Central_path(f,A,b,tolerance,fig);
disp("########### Central_path fixed results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Central_path_fixed(f,A,b,tolerance,fig);
disp("########### Central_path fixed (better compution effciency func results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Central_path_adaptive(f,A,b,tolerance,fig);
disp("########### Central_path adaptive results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=central_path_adaptive_2(f,A,b,tolerance,fig);
disp("########### Central_path adaptive2 better compution func results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

%% %Test 2 
clear
clc
close all
f=[-30,-20];
A=[2,1;1,3];
b=[8,8];
n=.91;
tolerance=.0001;
fig=1;             %%flag to display figures 
[X_values,Fmin_values,X_S]=Mehrotra_Predictor(f,A,b,n,tolerance,fig);
disp("########### Moherete predictor results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Mehrotra_Predictor_2(f,A,b,n,tolerance,fig);
disp("########### Moherete 2 (better compution effciency predictor results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))


[X_values,Fmin_values,X_S]=Central_path(f,A,b,tolerance,fig);
disp("########### Central_path fixed results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Central_path_fixed(f,A,b,tolerance,fig);
disp("########### Central_path fixed (better compution effciency func results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Central_path_adaptive(f,A,b,tolerance,fig);
disp("########### Central_path adaptive results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=central_path_adaptive_2(f,A,b,tolerance,fig);
disp("########### Central_path adaptive2 better compution func results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

%% Test 3
clear
clc
close all
f=[-3,-2,-5];
A=[1 2 1;3 0 2;1 4 0];
b=[430 460 420];

n=.91;
tolerance=.0001;
fig=1;             %%flag to display figures 
[X_values,Fmin_values,X_S]=Mehrotra_Predictor(f,A,b,n,tolerance,fig);
disp("########### Moherete predictor results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Mehrotra_Predictor_2(f,A,b,n,tolerance,fig);
disp("########### Moherete 2 (better compution effciency predictor results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))


[X_values,Fmin_values,X_S]=Central_path(f,A,b,tolerance,fig);
disp("########### Central_path fixed results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Central_path_fixed(f,A,b,tolerance,fig);
disp("########### Central_path fixed (better compution effciency func results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=Central_path_adaptive(f,A,b,tolerance,fig);
disp("########### Central_path adaptive results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))

[X_values,Fmin_values,X_S]=central_path_adaptive_2(f,A,b,tolerance,fig);
disp("########### Central_path adaptive2 better compution func results #####    ")
disp("X_values");
disp(X_values(end,:))
disp("F Min value ")
disp(Fmin_values(end))





















