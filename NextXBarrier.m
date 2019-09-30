%% Definition of the function

function xkp1 = NextX(xk, u, Size)
    xkp1 = zeros(1, length(u));
    for i = 1 : length(u)
        %nextCandidate =  round(2*sqrt(abs(xk - u(i))) + 2*xk + u(i), 0);
        %this is a dummy function to try
        nextCandidate = round(xk + u(i), 0);

        if nextCandidate > Size

            xkp1(i) = Size;
            
        elseif nextCandidate < 1
            
            xkp1(i) = 1;

        else
            xkp1(i) = nextCandidate;
        end
        
        if xk<40 && xkp1>=40
            xkp1=39;
        end
        if xk>40 && xkp1<=40
            xkp1=41;
        end

    end
end

