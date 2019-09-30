%Antoine dice che loro comunque discretizzano l'input. Quindi lo faccio
%anche io per trovare il minimo.
%It seems to be working, even when I introduce the barriers in the dynamics
%so that some state cannot reach the safe state.

clear variables
close all

Size  = 100;
PLOT = false;

%For robustness
lambda = 0.8;

Iterations = 100;

U = -3:3;

%time for the evolution of the dynamics
t = 1:100;

%% Initial condition (height)

h = zeros(Size, 1);

for i = 1 : Size
%         
%         %linear h
%         h(i, 1) = abs(-i + Size/2)-Size/5;
%         %sin h
        %with this function for h, when we introduce uncertainties, we see
        %that problems arise
        h(i, 1) = 100*sin(i/5);      
            % %Some cannot reach save zone
            % %Remember to use the barrier NEXT
            % h(i, 1) = 2+100/((i-20)^2-400);  

end


%% Update Vk (with discretized input)

Vk = zeros (Size, Iterations);

PossibleVkp1 = zeros(length(U), 1);

for i = 1:Size
    Vk(i, 1)=h(i);
end

for iter = 2:Iterations
    
    for xk = 1:Size
    %calculate the possible vk varying u
        for i = 1 : length(U)
           u = U(i);
           PossibleVkp1(i) = Vk(NextX(xk, u, Size), iter-1);

        end

    %choose the minimum one
    Vk(xk, iter) = max ( h(xk), min(PossibleVkp1)/lambda);
    
    end

end

%% Control and dynamics 

Xdyn = zeros (Size, length(t));
Hdyn = zeros (Size, length(t));
PossibleVstar = zeros(length(U), 1);
Vstar = Vk(:, Iterations);

if PLOT
    plot(1:Size, h)
    hold on
end

for i = 1:Size
    Xdyn(i, 1)=i;
    Hdyn(i, 1)=h(i);
end


for iter = 2:length(t)
    
    for xk = 1:Size
         for i = 1 : length(U)
           u = U(i);
           PossibleVstar(i) = Vstar(NextX(Xdyn(xk, iter-1), u, Size));
           
         end
    [~ ,uindex]= min(PossibleVstar);
%     %This is for the discrete normal case
%     Xdyn(xk, iter) = NextX(Xdyn(xk, iter-1), U(uindex), Size);
    %This is for the noisy (or uncertain) case
    Xdyn(xk, iter) = NextUncertain(Xdyn(xk, iter-1), U(uindex), Size);
    Hdyn (xk, iter) = h(Xdyn(xk, iter));
    
    
    end
    
if PLOT
    plot(Xdyn(10, iter-1),Hdyn (10, iter-1) , 'ro');
    
end
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




