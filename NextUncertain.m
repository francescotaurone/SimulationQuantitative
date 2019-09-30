%% Definition of the function, this is in order to create the noise
%We deal with it using the robust lambda method

function xkp1 = NextUncertain(xk, u, Size)
    xkp1 = zeros(1, length(u));
    for i = 1 : length(u)
        
        %Round so to get the index
        nextCandidate = round(xk + 2 + u(i), 0); 

        if nextCandidate > Size

            xkp1(i) = Size;
            
        elseif nextCandidate < 1
            
            xkp1(i) = 1;

        else
            xkp1(i) = nextCandidate;
        end

    end
end

