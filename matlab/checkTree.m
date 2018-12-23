function checkTree()
global counter;
counter = 0;
A = load('tree.csv');

recFun(A , 0 , 0 , 2 , A(1,:) , 2);
counter
for i = 2:length(A(:,1))
    dist = (A(1,1)-A(i,1))^2 + (A(1,2) - A(i,2))^2 ;
end
end

function c = check(vp , p , median)

   dist = 0;
   c = 2;
   for i = 1:length(vp)
       dist = dist + (vp(i) - p(i))^2;
   end
   [dist , median];
   if dist  < median
       c = 0;
   elseif dist > median
       c = 1;
   end


end
function recFun(A , i , level , maxLevel , vpmedian ,lowhighroot) %012 low high root
    global counter;
    counter = counter + 2;
    if(level+1 < maxLevel)
        if lowhighroot >-1
            'roooot';
            A(i+1,:);
        end
       indexlow = (i+1) + 1;
       indexhigh = (i+1) + 2^(maxLevel - level);
       c1 = check(vpmedian(1:2),A(indexlow,1:2) , vpmedian(3));
       c2 = check(vpmedian(1:2),A(indexhigh,1:2) , vpmedian(3));
       if lowhighroot == 0
           if c1 == 1
               'wrong'
           end
           if c2 == 1
               'wrong'
           end
       elseif lowhighroot == 1
           if c1 == 0
               'wrong'
           end
           if c2 == 0
               'wrong'
           end
       else
           if c1 == 1
               'wrong'
           end
           if c2 == 0
               'wrong'
           end
       end
   
       recFun(A , i+1 , level+1 , maxLevel , vpmedian , 0);
       recFun(A , i+1+2^(maxLevel - level) , level+1 , maxLevel , vpmedian , 1);
        
       
       recFun(A , i+1 , level+1 , maxLevel , A(i+2,:) , 2);
       recFun(A , i+1+2^(maxLevel - level) , level+1 , maxLevel , A(i+1+2^(maxLevel - level),:) , 2);
        
        
    end
    




end