function out = circle_inters()
    %close all
    point_color = 'k';
    nbr_color = 'k';
    th = 0:pi/500:2*pi;
    vpIndex = 1;
    A = [1 1 1; 2 1 0.75; 3 2 1.75; 1 3 1.5];
    A = load('tree.csv');
    maxLevel = log2(length(A(:,1)) + 1);
    recRemove(A , 0 , [] , [] , 0, maxLevel-1,th);
    plot(A(:,1) , A(:,2) , '.' ,'MarkerSize' , 10 , 'MarkerFaceColor' , point_color ,'MarkerEdgeColor',point_color)
    nbrs = vpKnn(A(vpIndex,:));
    hold on;
    plot(nbrs(:,1) , nbrs(:,2) , '+' ,'MarkerSize' , 10 , 'MarkerFaceColor' , nbr_color ,'MarkerEdgeColor',nbr_color);
    hold on;
    plot(A(vpIndex,1) , A(vpIndex,2) , 'o' ,'MarkerSize' , 10 , 'MarkerFaceColor' , point_color ,'MarkerEdgeColor',point_color);

end


function [x,y] = circlePoints(center , radius , th)
    x = center(1) + radius*cos(th);
    y = center(2) + radius*sin(th);
end

function [x,y] = removeIntersectingPointsIns(x, y , center , radius)
    
    x1 = [];
    y1 = [];
    for i = 1:length(x)
        if i < length(x)
            if (x(i)-center(1))^2 + (y(i)-center(2))^2 < radius^2
                'removed';
                x1 = [x1 , x(i)];
                y1 = [y1 , y(i)];
            end
        end
    end
   x = x1;
   y = y1;

end
function [x,y] = removeIntersectingPointsOuts(x, y , center , radius)
    x1 = [];
    y1 = [];
    for i = 1:length(x)
        if i < length(x)
            if (x(i)-center(1))^2 + (y(i)-center(2))^2 > radius^2
                'removed';
                x1 = [x1 , x(i)];
                y1 = [y1 , y(i)];
            end
        end
    end
   x = x1;
   y = y1;

end
function  x = recRemove(A , i  , ins , outs , level , maxLevel , th)
    i;
    level;
    Colors = [[0.5,0,1] ; [1 , 0.2 , 0.8] ; [1 , 0.4 , 0.6] ; [0.5 , 0.6 , 0.4] ; [0 , 0.3 , 0.7]; [0.4 , 0.3 , 0.8]; [1 , 0.2 , 0.8] ; [1 , 0.4 , 0.6] ; [0.5 , 0.6 , 0.4] ; [0 , 0.3 , 0.7]; [0.4 , 0.3 , 0.8]];
    Colors(level+1);
    if maxLevel > level
        [x,y] = circlePoints([A(i+1,1) , A(i+1,2)] , A(i+1,3) , th);
        if length(outs)>0
         for j = 1:length(outs(:,1))
             "outs";
             j;
             outs;
            [x,y] = removeIntersectingPointsOuts(x,y,[outs(j,1) , outs(j,2)] , outs(j,3));
         end
        end
        if length(ins) > 0
          for j = 1:length(ins(:,1))
              ins;
              j;
            [x,y] = removeIntersectingPointsIns(x,y,[ins(j,1) , ins(j,2)] , ins(j,3));
          end
        end
          plot(x , y , '.','MarkerEdgeColor',Colors(level+1,:),'MarkerFaceColor',Colors(level+1,:));
          hold on;
         recRemove(A , i + 1 , [ins ; A(i+1,:)] , outs , level+1 , maxLevel,th); 
         recRemove(A , i + 2^(maxLevel - level) , ins , [outs ; A(i+1,:)] , level+1 , maxLevel,th);
        
        
    end


end