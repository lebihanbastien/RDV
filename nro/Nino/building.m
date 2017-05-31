function Lam = building(choice, sol, np, it_max)

for i = 1:np-1
     
     if choice(1)
        it_CW  = sol(i).CW.it;    % to check convergence it is suff to check these variables
        Lam_CW = sol(i).CW.Lam(:,1:12);
     end
     
     if choice(2)
        it_SL  = sol(i).SL.it;
        Lam_SL = sol(i).SL.Lam(:,1:12);
     end
     
     if choice(3)
        it_LR  = sol(i).LR.it;
        Lam_LR = sol(i).LR.Lam(:,1:12);
     end
     
     if choice(4) 
        it_cont  = sol(i).cont.it;
        Lam_cont = sol(i).cont.Lam(:,1:12);
     end
     
%% Clohessy - Wiltshire

     if choice(1) && ~choice(2)
        if it_CW < it_max      
           Lam_temp = Lam_CW;
        else
           error('C-W did not converge for transfer %d',i)
        end
     end  
           
%% Straight Line     
     
       if choice(2) && ~choice(1)
        if it_SL < it_max      %convergence
           Lam_temp = Lam_SL;
        else
           error('SL did not converge for transfer %d',i)
        end
       end  

%% Luquette       
     if choice(3) && ~choice(1)
        if it_LR < it_max      %convergence
           Lam_temp = Lam_LR;
        else
           error('Luquette did not converge for transfer %d',i)
        end
     end  
           
           
%% Multiple Choice

     if choice(1) && choice(2) && choice(3) 
        if it_LR < it_max      %convergence
           sprintf('Lam for Transfer %d from Luquette',i)
           Lam_temp = Lam_LR;
           
        elseif it_CW < it_max       % no conv + multiple choice
               sprintf('Lam for Transfer %d from C-W',i)
               Lam_temp = Lam_CW;
               
        else 
               if it_SL < it_max
                  sprintf('Lam for Transfer %d from SL',i)
                  Lam_temp = Lam_SL;
               else
                   sprintf('All methods failed to converge for transfer %d',i)
                   disp('Continuation Algorithm')
                   error('puzza')
               end 
        end
     end
    
     
 %% Continuation    
     if choice(4) 
        if it_cont < it_max      %convergence
           Lam_temp = Lam_cont;
        else
           error('Continuation did not converge for transfer %d',i)
        end
     end       
     
 
 %% Final output
 
 Lam{i} = Lam_temp; 
 
 end
 