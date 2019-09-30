%Antoine dice che loro comunque discretizzano l'input. Quindi lo faccio
%anche io per trovare il minimo.
%It seems to be working, even when I introduce the barriers in the dynamics
%so that some state cannot reach the safe state.

clear variables
close all

Size  = 10;
PLOT = false;

%For robustness
lambda = 0.8;

Iterations = 50;

U1 = -5:5;
U2 = -5:5;
%time for the evolution of the dynamics
t = 1:100;

%% Initial condition (height)

h = zeros(Size, Size, 1);

for i = 1 : Size
for j = 1 : Size
%         
%         %linear h
%         h(i,j, 1) = abs(-(i+j) + Size/2)-Size/5;
%         %sin h 
%         %with this function for h, when we introduce uncertainties, we see
%         %that problems arise
%         h(i,j, 1) = 100*sin((i+j)/5);      
            % %Some cannot reach save zone
            % %Remember to use the barrier NEXT
            % h(i, 1) = 2+100/((i-20)^2-400);
            %Quadratic one
              
        h(i,j, 1) = i^2 + j^2-100;
end
end


%% Update Vk (with discretized input)

Vk = zeros (Size,Size, Iterations);

PossibleVkp1 = zeros(length(U1),length(U2), 1);

for i = 1:Size
    for j = 1:Size
    Vk(i,j, 1)=h(i,j);
    end
end
for iter = 2:Iterations
    
    for x1k = 1:Size
    for x2k = 1:Size
    %calculate the possible vk varying u
        for i = 1 : length(U1)
        for j = 1 : length(U2)
           u1 = U1(i);
           u2 = U2(j);
           [x1kp1, x2kp1] = NextX2D(x1k, x2k, u1, u2, Size);
           PossibleVkp1(i, j) = Vk(x1kp1, x2kp1, iter-1);

        end
        end
    %choose the minimum one
    Vk(x1k,x2k, iter) = max ( h(x1k, x2k), min(min(PossibleVkp1(:,:)))/lambda);
    
    end
    end
end

%% Control and dynamics 

X1dyn = zeros (Size,Size, length(t));
X2dyn = zeros (Size,Size, length(t));
Hdyn = zeros (Size,Size, length(t));
PossibleVstar = zeros(length(U1),length(U2), 1);
Vstar = Vk(:,:, Iterations);

% if PLOT
%     plot(1:Size, h)
%     hold on
% end

for i = 1:Size
for j = 1:Size
    X1dyn(i,j, 1)=i;
    X2dyn(i,j, 1)=j;
    Hdyn(i,j, 1)=h(i,j);
end
end

for iter = 2:length(t)
    
    for x1k = 1:Size
    for x2k = 1:Size
         for i = 1 : length(U1)
         for j = 1 : length(U2)
           u1 = U1(i);
           u2 = U2(i);
           [x1kp1, x2kp1] = NextX2D(x1k, x2k, u1, u2, Size);
           PossibleVstar(i,j) = Vstar(x1kp1, x2kp1);
           
         end
         end
   
    [u1index, u2index]=find(PossibleVstar == min(min(PossibleVstar(:,:))));
    
    %This is for the discrete normal case
    pastx1 =X1dyn(x1k,x2k, iter-1);
    pastx2 =X2dyn(x1k,x2k, iter-1);
     %WATCHOUT! HERE i ADDED THE 1S JUST IN CASE OF MULLTIPLE APPEARANCES
    [X1dyn(x1k,x2k, iter), X2dyn(x1k,x2k, iter)] = NextX2D(pastx1,pastx2, U1(u1index(1)),U2(u2index(1)), Size);
%     %This is for the noisy (or uncertain) case
%     Xdyn(x1k, iter) = NextUncertain(Xdyn(x1k, iter-1), U(uindex), Size);
    newx1 =X1dyn(x1k,x2k, iter);
    newx2 =X2dyn(x1k,x2k, iter);
    Hdyn (x1k,x2k, iter) = h(newx1, newx2);
    
    
    end
    end
    
% if PLOT
%     plot(Xdyn(10, iter-1),Hdyn (10, iter-1) , 'ro');
%     
% end
end

%% Update Vk (with no discretization, still not working)

%This is the one using fmincon, that is not the best iun this context apparently 
% xk = 1;
% 
% 
% f = @(u) FunctionToMin(xk, u, Size, h);
% %f = @(u) DummyFunction(u);
% u0 = 110;
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = -InputConstr;
% ub = InputConstr;
% usolution = fmincon(f,u0,A,b,Aeq,beq,lb,ub)




