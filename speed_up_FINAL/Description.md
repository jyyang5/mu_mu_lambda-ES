### speed_up_FINAL



##### Analysis
T = [noGP,GP(N=10),GP(N=100),GP(\infty)];
speed_up_self = T(:,2:4)./repmat(T(:,1),1,3);
- Each row one strategy ((3/3,10),(5/5,20),(10/10,40),(1+1))

##### Surrogate assisted mml-ES CSA
T_1 = [1+1,(3/3,10),(5/5,20),(10/10,40)]; 
speed_up_1 = repmat(T_1(:,1),1,3)./T_1(:,2:4);

- Each row one test function (linear,quadratic,cubic,Schwefel,quartic)

##### Surrogate assisted mml-ES CSA(plus-selection)
T_2 = [1+1,(3/3,10),(5/5,20),(10/10,40)];  
speed_up_2 = repmat(T_2(:,1),1,3)./T_2(:,2:4);
- Each row one test function (linear,quadratic,cubic,Schwefel,quartic)
