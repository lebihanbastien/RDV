 %Graphical comparison of final Lambert's solution to verify uniqueness
 figure
 for i=1:np-1
 p1 = plot3(Transfer_CW(i).it(end).yv(:,1),Transfer_CW(i).it(end).yv(:,2),Transfer_CW(i).it(end).yv(:,3),'r');
 hold on
 grid on
 p2 = plot3(Transfer_SL(i).it(end).yv(:,1),Transfer_SL(i).it(end).yv(:,2),Transfer_SL(i).it(end).yv(:,3),'go');
 p3 = plot3(Transfer_LR(i).it(end).yv(:,1),Transfer_LR(i).it(end).yv(:,2),Transfer_LR(i).it(end).yv(:,3),'ko');
 end


 title('Uniqueness of Lambert''s solution in the LVLH frame, LVLH dynamics')
 legend('CW','SL','LR')
 xlabel('x_{syn}') 
 ylabel('y_{syn}')
 zlabel('z_{syn}')


%first guesses and Lambert solution synodic 
figure
 for i=1:np-1
 p1 = plot3(Transfer_CW(i).it(1).yv(:,1),Transfer_CW(i).it(1).yv(:,2),Transfer_CW(i).it(1).yv(:,3),'r');
 hold on
 grid on
 %plot3(Transfer_CW(i).it(end).yv(:,1),Transfer_CW(i).it(end).yv(:,2),Transfer_CW(i).it(end).yv(:,3),'g')
 p2 = plot3(Transfer_SL(i).it(1).yv(:,1),Transfer_SL(i).it(1).yv(:,2),Transfer_SL(i).it(1).yv(:,3),'m');
 %plot3(Transfer_SL(i).it(end).yv(:,1),Transfer_SL(i).it(end).yv(:,2),Transfer_SL(i).it(end).yv(:,3),'g');
 p3 = plot3(Transfer_LR(i).it(1).yv(:,1),Transfer_LR(i).it(1).yv(:,2),Transfer_LR(i).it(1).yv(:,3),'k');
 p4 = plot3(Lam_LVLH_hyb(i).transfer.yv(:,1),Lam_LVLH_hyb(i).transfer.yv(:,2),Lam_LVLH_hyb(i).transfer.yv(:,3),'g');
 p6(i) = plot3(Lam_LVLH_hyb(i).transfer.yv(end,1),Lam_LVLH_hyb(i).transfer.yv(end,2),Lam_LVLH_hyb(i).transfer.yv(end,3),'bo');

 end
 p5 = plot3(r_f(1,1), r_f(1,2), r_f(1,3),'co');
 p7 = plot3(0,0,0,'ko');
 
 if np>2  
   for i=1:np-2
      str(i) = cellstr(sprintf('%s %d','HP',i));
   end 
   pp = [p1 p2 p3 p4 p5 p6(1:np-2) p7];
   nom = ['C-W' 'SL' 'LR' 'Lambert opt' 'Start' str 'Docking'];
   legend(pp,nom)
else
   legend([p1 p2 p3 p4 p5 p7],'C-W', 'SL', 'LR', 'Lambert opt', 'Start', 'Docking')
 end

 title('First Guesses and Lambert Solution in the LVLH Frame with LVLH dynamics')
 xlabel('x_{LVLH}') 
 ylabel('y_{LVLH}')
 zlabel('z_{LVLH}')
 
% %Error
% figure
%  for i=1:np-1
%  plot(i,err_CW(i).tr(1).*cr3bp.L,'ro')
%  hold on
%  plot(i,err_SL(i).tr(1).*cr3bp.L,'bo')
%  plot(i,err_LR(i).tr(1).*cr3bp.L,'go')
%  end
%  
%  legend('CW','SL','LR')
%  title('Initial Dimensional Error for each Transfer with LVLH Dynamics')
%  
%  xlabel('Transfer')
%  ylabel('Error [km]')

 
%  Lambert solution in LVLH
%  figure
%  for i=1:np-1
%      plot3(Lam_LVLH(i).yv(:,1),Lam_LVLH(i).yv(:,2),Lam_LVLH(i).yv(:,3))
%      hold on
%      grid on
%  end
%   
%  view(180,180)
%  plot3(0,0,0,'ro')
%  
%  title('Lambert Arc in LVLH')
%  xlabel('x')
%  ylabel('y')
%  zlabel('z')
 
 
%  %RDV in synodic frame
%  figure
%  plot3(nro_T.yv(:,1),nro_T.yv(:,2),nro_T.yv(:,3),'b')  %Target 
%  hold on
%  axis equal
%  plot3(nro_C.yv(:,1),nro_C.yv(:,2),nro_C.yv(:,3),'g')  %Chaser 
%  plot3(x0_T(1), x0_T(2), x0_T(3),'ro')
%  plot3(x0_C(1), x0_C(2), x0_C(3),'ko') 
%  for i=1:np-1
%  plot3(Transfer_LR(i).it(end).yv(:,7)+Transfer_LR(i).it(end).yv(:,1), Transfer_LR(i).it(end).yv(:,8)+Transfer_LR(i).it(end).yv(:,2), Transfer_LR(i).it(end).yv(:,9)+Transfer_LR(i).it(end).yv(:,3),'r')
%  plot3(Transfer_LR(i).it(end).yv(end,7)+Transfer_LR(i).it(end).yv(end,1), Transfer_LR(i).it(end).yv(end,8)+Transfer_LR(i).it(end).yv(end,2), Transfer_LR(i).it(end).yv(end,9)+Transfer_LR(i).it(end).yv(end,3),'ko')
%  end
% R_Moon = 1737/cr3bp.L; 
% MOON = imread('Moon.jpg','jpg');
% props.FaceColor='texture';
% props.EdgeColor='none';
% props.FaceLighting='phong';
% props.Cdata = MOON;
% Center = [1 - cr3bp.mu; 0; 0];
% [XX, YY, ZZ] = ellipsoid(-Center(1),Center(2),Center(3),R_Moon,R_Moon,R_Moon,30);
% surface(-XX, -YY, -ZZ,props);
% 
%  %legend('Target Orbit','Chaser Orbit','Target at t0','Chaser at t0','Corrected Trajectory','Docking')
%  title('Nino'' s RDV (Synodic Frame)')