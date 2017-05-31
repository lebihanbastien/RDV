function [] = plot_3_choice(choice, sol, Lam1, r_f_syn, np, cr3bp, nro_T, nro_C, x0_T, x0_C)

%% plots for choice = [1 1 1 0]

%Graphical comparison of final Lambert's solution to verify uniqueness
  if choice(1) && choice(2)  && choice(3) 
%  figure     
%  for i=1:np-1
%      Lam_CW = sol(i).CW.Lam;
%      Lam_SL = sol(i).SL.Lam;
%      Lam_LR = sol(i).LR.Lam;
%      
%      p1 = plot3(Lam_CW(:,1),Lam_CW(:,2),Lam_CW(:,3),'r');
%      hold on
%      grid on
%      p2 = plot3(Lam_SL(:,1), Lam_SL(:,2), Lam_SL(:,3),'gd');
%      p3 = plot3(Lam_LR(:,1), Lam_LR(:,2), Lam_LR(:,3),'ko');
%  end
% 
% 
%  title('Uniqueness of solution in the Target-Centered Synodic Frame')
%  legend('CW','SL','L')
%  xlabel('x_{syn}') 
%  ylabel('y_{syn}')
%  zlabel('z_{syn}')
%  
%  %first guesses and Lambert solution synodic 
% figure
%  for i=1:np-1
%      fg_CW = sol(i).CW.fg;
%      fg_SL = sol(i).SL.fg;
%      fg_LR = sol(i).LR.fg;
%      Lam   = Lam1{i};
% 
% 
%      p1 = plot3(fg_CW(:,1), fg_CW(:,2), fg_CW(:,3),'r');
%      hold on
%      grid on
%      p2 = plot3(fg_SL(:,1), fg_SL(:,2), fg_SL(:,3),'m');
%      p3 = plot3(fg_LR(:,1), fg_LR(:,2), fg_LR(:,3),'k');
%      p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
%      p6(i) = plot3(Lam(end,1), Lam(end,2), Lam(end,3),'bo');
% 
%  end
%  p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
%  p7 = plot3(0,0,0,'ko');
%    
% if np>2  
%    for i=1:np-2
%       str(i) = cellstr(sprintf('%s %d','HP',i));
%    end 
%    pp = [p1 p2 p3 p4 p5 p6(1:np-2) p7];
%    nom = ['C-W' 'SL' 'L' 'Lambert opt' 'Start' str 'Docking'];
%    legend(pp,nom)
% else
%    legend([p1 p2 p3 p4 p5 p7],'C-W', 'SL', 'L', 'Lambert', 'Start', 'Docking')
% end
%  
%  title('First Guesses and Lambert Solution, Synodic Frame')
%  xlabel('x_{syn}') 
%  ylabel('y_{syn}')
%  zlabel('z_{syn}')
%  
 
%Dimensional Error
 figure
 for i=1:np-1
     err_CW = sol{i}.CW.err1;
     err_SL = sol{i}.SL.err1;
     err_LR = sol{i}.LR.err1;

     semilogy(i,err_CW*cr3bp.L,'ro')
     hold on
     grid on
     semilogy(i,err_SL*cr3bp.L,'mo')
     semilogy(i,err_LR*cr3bp.L,'ko')
 end

 set(gca,'XTick',[1:i] );
 legend('CW_{syn}','SL_{syn}','L_{syn}')
 title('Init Dim Error for each Transfer, Syn Dynamics')
 
 xlabel('Transfer')
 ylabel('Error [km]')
 
 end
 
%% plots for choice = [1 0 0 0] (C-W)
 
if choice(1) && ~choice(2) 
%first guesses and Lambert solution synodic 
   figure
   for i = 1:np-1
       fg_CW = sol{i}.CW.fg;
       Lam   = Lam1{i};

       p1 = plot3(fg_CW(:,1), fg_CW(:,2), fg_CW(:,3),'r');
       hold on
       grid on
       p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
       p6(i) = plot3(Lam(end,1), Lam(end,2), Lam(end,3),'bo');

   end
   p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
   p7 = plot3(0,0,0,'ko');
   
   if np>2  
      for i=1:np-2
          str(i) = cellstr(sprintf('%s %d','HP',i));
      end 
      pp = [p1 p4 p5 p6(1:np-2) p7];
      nom = ['C-W' 'Lambert opt' 'Start' str 'Docking'];
      legend(pp,nom)
   else
      legend([p1 p4 p5 p7],'C-W', 'Lambert', 'Start', 'Docking')
   end
 
   title('CW and Lambert Solution, Synodic Frame')
   xlabel('x_{syn}') 
   ylabel('y_{syn}')
   zlabel('z_{syn}')
 
 
 %Dimensional Error
 figure
 for i=1:np-1
     err_CW = sol{i}.CW.err(1);
     semilogy(i, err_CW*cr3bp.L,'ro')
     hold on
     grid on
 end

 set(gca,'XTick',[1:i] );
 legend('CW_{syn}')
 title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
 xlabel('Transfer')
 ylabel('Error [km]')
 
end

%% plots for choice = [0 1 0 0] (SL)


if choice(2)  && ~choice(1) 
%first guesses and Lambert solution synodic 
figure
 for i=1:np-1
     fg_SL = sol{i}.SL.fg;
     Lam   = Lam1{i};
       
     p2 = plot3(fg_SL(:,1), fg_SL(:,2), fg_SL(:,3),'m');
     hold on
     grid on
     p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
     p6(i) = plot3(Lam(end,1), Lam(end,2), Lam(end,3),'bo');

 end
 p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
 p7 = plot3(0,0,0,'ko');
   
if np>2  
   for i=1:np-2
      str(i) = cellstr(sprintf('%s %d','HP',i));
   end 
   pp = [p2 p4 p5 p6(1:np-2) p7];
   nom = ['SL' 'Lambert opt' 'Start' str 'Docking'];
   legend(pp,nom)
else
   legend([ p2 p4 p5 p7],'SL', 'Lambert', 'Start', 'Docking')
end
 
 title('SL and Lambert Solution, Synodic Frame')
 xlabel('x_{syn}') 
 ylabel('y_{syn}')
 zlabel('z_{syn}')
 
 
 %Dimensional Error
 figure
 for i=1:np-1
     err_SL = sol{i}.SL.err(1);
     semilogy(i,err_SL*cr3bp.L,'mo')
     hold on
     grid on
 end

 set(gca,'XTick',[1:i] );
 legend('SL_{syn}')
 title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
 xlabel('Transfer')
 ylabel('Error [km]')
 
end
 

%% plots for choice = [0 0 1 0] (Luquette)

if choice(3)  && ~choice(1) 
%first guesses and Lambert solution synodic 
figure
 for i=1:np-1
     fg_LR = sol{i}.LR.fg;
     Lam   = Lam1{i};
     
     p3 = plot3(fg_LR(:,1), fg_LR(:,2), fg_LR(:,3),'k');
     hold on
     grid on
     p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
     p6(i) = plot3(Lam(end,1), Lam(end,2), Lam(end,3),'bo');

 end
 p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
 p7 = plot3(0,0,0,'ko');
   
if np>2  
   for i=1:np-2
      str(i) = cellstr(sprintf('%s %d','HP',i));
   end 
   pp = [p3 p4 p5 p6(1:np-2) p7];
   nom = ['L' 'Lambert opt' 'Start' str 'Docking'];
   legend(pp,nom)
else
   legend([p3 p4 p5 p7],'L', 'Lambert', 'Start', 'Docking')
end
 
 title('Luquette and Lambert Solution, Synodic Frame')
 xlabel('x_{syn}') 
 ylabel('y_{syn}')
 zlabel('z_{syn}')
 
  %Dimensional Error
 figure
 for i=1:np-1
     err_LR = sol{i}.LR.err(1);
     semilogy(i,err_LR*cr3bp.L,'mo')
     hold on
     grid on
 end

 set(gca,'XTick',[1:i] );
 legend('L_{syn}')
 title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
 xlabel('Transfer')
 ylabel('Error [km]')
 
end
 

%% plots for choice = [0 0 0 1] (Continuation)

if choice(4) 
%first guesses and Lambert solution synodic 
figure
 for i=1:np-1
     fg_cont = sol{i}.cont.fg;
     Lam   = Lam1{i};
     
     p3 = plot3(fg_cont(:,1), fg_cont(:,2), fg_cont(:,3),'k');
     hold on
     grid on
     p4 = plot3(Lam(:,1), Lam(:,2), Lam(:,3),'g');
     p6(i) = plot3(Lam(end,1), Lam(end,2), Lam(end,3),'bo');

 end
 p5 = plot3(r_f_syn(1,1), r_f_syn(1,2), r_f_syn(1,3),'co');
 p7 = plot3(0,0,0,'ko');
   
if np>2  
   for i=1:np-2
      str(i) = cellstr(sprintf('%s %d','HP',i));
   end 
   pp = [p3 p4 p5 p6(1:np-2) p7];
   nom = ['Continuation' 'Lambert opt' 'Start' str 'Docking'];
   legend(pp,nom)
else
   legend([p3 p4 p5 p7],'L', 'Lambert', 'Start', 'Docking')
end
 
 title('Continuation and Lambert Solution, Synodic Frame')
 xlabel('x_{syn}') 
 ylabel('y_{syn}')
 zlabel('z_{syn}')
 
  %Dimensional Error
 figure
 for i=1:np-1
     err_cont = sol{i}.cont.err(1);
     semilogy(i, err_cont*cr3bp.L,'mo')
     hold on
     grid on
 end

 set(gca,'XTick',[1:i] );
 legend('Continuation_{syn}')
 title('Initial Dimensional Error for each Transfer, Synodic Dynamics')
 
 xlabel('Transfer')
 ylabel('Error [km]')
 
end

%% plots plotted in any case

% RDV in synodic frame
figure
p1 = plot3(nro_T.yv(:,1),nro_T.yv(:,2),nro_T.yv(:,3),'g');  %Target 
hold on
axis equal
grid on
p2 = plot3(nro_C.yv(:,1),nro_C.yv(:,2),nro_C.yv(:,3),'k');  %Chaser 
p3 = plot3(x0_T(1), x0_T(2), x0_T(3),'mo');
p4 = plot3(x0_C(1), x0_C(2), x0_C(3),'co','LineWidth',1.5); 
 for i=1:np-1
     Lam = Lam1{i};
     p5 = plot3(Lam(:,7) + Lam(:,1), Lam(:,8) + Lam(:,2), Lam(:,9) + Lam(:,3),'r','LineWidth',1.5);
     p6(i) = plot3(Lam(end,7) + Lam(end,1), Lam(end,8) + Lam(end,2), Lam(end,9) + Lam(end,3),'bo','LineWidth',1.5);
 end
p7 = plot3(Lam(end,7) + Lam(end,1), Lam(end,8) + Lam(end,2), Lam(end,9) + Lam(end,3),'ko','LineWidth',1.5);

R_Moon = 1737/cr3bp.L; 
MOON = imread('moon.jpg','jpg');
props.FaceColor='texture';
props.EdgeColor='none';
props.FaceLighting='phong';
props.Cdata = MOON;
Center = [1 - cr3bp.mu; 0; 0];
[XX, YY, ZZ] = ellipsoid(-Center(1),Center(2),Center(3),R_Moon,R_Moon,R_Moon,30);
surface(-XX, -YY, -ZZ,props);

 if np>2  
   for i=1:np-2
      str(i) = cellstr(sprintf('%s %d','HP',i));
   end 
   pp = [p1 p2 p3 p4 p5 p6(1:np-2) p7];
   nom = ['Target Orbit' 'Chaser Orbit' 'Target initial position' 'Chaser initial position' 'RDV trajectory' str 'Docking'];
   legend(pp,nom)
else
   legend([p1 p2 p3 p4 p5 p7],'Target Orbit', 'Chaser Orbit', 'Target initial position', 'Chaser initial position', 'RDV trajectory' ,'Docking')
  end

 title('Rendezvous in the Synodic Frame')
 xlabel('x_{syn}')
 ylabel('y_{syn}')
 zlabel('z_{syn}')

