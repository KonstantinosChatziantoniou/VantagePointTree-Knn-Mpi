function nbrsArr = vpKnn(vpIn)
%% globals and init
global tau;
global nbrs;
global nbrsLen;
global k;       %const
global vp;      %const
global maxLevel;%const
global A;       %const(cols 1 , 2)  not const col3

tau = Inf; %no points in nbrs yet
nbrs = [];
nbrsLen = 0;
A = load('tree.csv');
maxLevel = log2(length(A(:,1))+1)-1;
vp = vpIn;
k = 3;




%% searching
recKnn(0 , 0);

%% out
nbrsArr = nbrs;

vpKdists = (nbrs(:,1) - vp(1)).^2 + (nbrs(:,2) - vp(2)).^2;
vpKdists = sqrt(vpKdists);
vpKdists = sort(vpKdists)
%% bruteforce KNN
nbrsArr2 = [];
len2 = 1;
maxNb = Inf;
for i = 1:length(A(:,1))
    A(i,3) = sqrt((A(i,1) - vp(1))^2 + (A(i,2) - vp(2))^2);
end

b = sort(A(:,3));
temp = 1;
if b(1) < 0.0000001
    temp = 2;
end
realMins = b(temp:(k+1))
if  realMins==vpKdists%abs(realMins-vpKdists(:,3))<0.0000001
    'yeeeeeeeeee'
else
    'fuuuuuuuuuu'
end




end
function recKnn(i , level)
    global tau;
    global maxLevel;
    global vp;
    global A;
    i;
    tau;
    level;
     %% if we are at leaf node check and return
    if level > maxLevel
        return
    end
    %% find distance from node
    dist = sqrt((A(i+1,1)-vp(1))^2 + (A(i+1,2) - vp(2))^2);
    
    [dist , tau]
    [dist , A(i+1,3)]
    %% if distance is lower than a threshold, the node is the same as the vp
    if dist < 0.0000001
        %do nothing
    %% compare distance from node  with tau
    else
        if dist < tau
            %if dist is lower than the median go to lower
            addNbr([A(i+1,1) , A(i+1,2) , dist]);
        end
    end
   
   
     %% searching next nodes
    if dist <= A(i+1 , 3) %if point is in the circle
        recKnn(i+1 , level+1);
        if dist + tau >= A(i+1 , 3)
            recKnn(i+2^(maxLevel - level) , level+1);
        end
    else
        recKnn(i+2^(maxLevel - level) , level+1);
        if dist - tau <= A(i+1,3)
            recKnn(i+1 , level+1);
        end
    end
            
        
end


function addNbr(point)
    global nbrsLen;
    global k;
    global nbrs;
    global tau;
    'asdasdasd'
    if nbrsLen < k
        nbrs = [nbrs ; point];
        nbrsLen = nbrsLen+1;
        'asdasdasd';
    else
        [m , i] = max(nbrs(:,3));
        if point(3) < m
            nbrs(i,:) = point;
        end
    end
    
    if nbrsLen == k
        tau = max(nbrs(:,3));
    end

end