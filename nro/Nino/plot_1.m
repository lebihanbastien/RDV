%Plotting
 figure
 for i=1:np-1
 p1 = plot3(state_CW(i).yv(:,1), state_CW(i).yv(:,2), state_CW(i).yv(:,3),'r');
 grid on
 hold on
 p2 = plot3(state_SL(i).yv(:,1), state_SL(i).yv(:,2), state_SL(i).yv(:,3),'g');
 p4(i) = plot3(state_CW(i).yv(end,1), state_CW(i).yv(end,2), state_CW(i).yv(end,3),'bo');
 end
 p3 = plot3(state_CW(1).yv(1,1), state_CW(1).yv(1,2), state_CW(1).yv(1,3),'co');
 p5 = plot3(0, 0, 0,'ko');
 
 
if np>2  
   for i=1:np-2
      str(i) = cellstr(sprintf('%s %d','HP',i));
   end 
   pp = [p1 p2 p3 p4(1:np-2) p5];
   nom = ['C-W' 'SL' 'Start' str 'Docking'];
   legend(pp,nom)
else
   legend([p1 p2 p3 p5],'C-W','SL','Start','Docking')
end
 
 title('CW and SL in the LVLH frame')
 xlabel('x_{LVLH}')
 ylabel('y_{LVLH}')
 zlabel('z_{LVLH}')
 
 
 
%  figure
%  for i=1:np-1
%  p1 = plot3(state_LR(i).yv(:,1), state_LR(i).yv(:,2), state_LR(i).yv(:,3),'r');
%  grid on
%  hold on
%  p3(i) = plot3(state_LR(i).yv(end,1), state_LR(i).yv(end,2), state_LR(i).yv(end,3),'bo');
% 
%  end
%  p2 = plot3(state_LR(1).yv(1,1), state_LR(1).yv(1,2), state_LR(1).yv(1,3),'co');
%  p4 = plot3(0, 0, 0,'ko');
%  
%  if np>2  
%    for i=1:np-2
%       str(i) = cellstr(sprintf('%s %d','HP',i));
%    end 
%    pp = [p1 p2 p3(1:np-2) p4];
%    nom = ['Luquette' 'Start' str 'Docking'];
%    legend(pp,nom)
% else
%    legend([p1 p2 p4],'Luquette','Start','Docking')
% end
% 
%  title('Luquette in the Target-Centered synodic frame')
%  xlabel('x_{syn}')
%  ylabel('y_{syn}')
%  zlabel('z_{syn}')