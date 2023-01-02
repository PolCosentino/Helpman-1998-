### Solve model
function  solveHLwCtyOpen_E(fund,dist,bord,bordc,nobs)


 #global alpha sigma LL LLwest LLeast;

 #xtic = tic();

 # Extract location characteristics from fundamentals matrix;
 # fund(:,1)=a; fund(:,2)=H; fund(:,3)=Iwest; fund(:,4)=Ieast;
 a = fund[:, 1]
 H = fund[:, 2]
 Iwest = fund[:, 3]
 Ieast = fund[:, 4]

 # convergence indicator;
 converge=0; 

 # Initialization based on a symmetric allocation;
 L_i = fill(LL/nobs, nobs, 1)
 w_i = fill(1.0, nobs, 1)


 # trade costs;
 dd = (dist .* bord .* bordc).^(1 - sigma)

 # ******************************************************;
 # **** START LOOP TO SOLVE FOR WAGES AND POPULATION ****;
 # ******************************************************;
 tradesh = 0
 dtradesh = 0

 x=1;
 while (x<200000)
     # Trade share;
     pwmat = (L_i .* (a .^ (sigma-1)) .* (w_i .^ (1-sigma))) * ones(1, nobs)
     nummat = dd .* pwmat
     denom = sum(nummat, dims=1)
     denommat = ones(nobs,1)*denom;
     tradesh = nummat ./ denommat

     # test;
     test=sum(tradesh, dims=1);
     mntest=mean(test);

     # Income equals expenditure equilibrium condition;
     income = w_i .* L_i
     expend=tradesh*income;    

     # domestic trade share;
     dtradesh = diag(tradesh)

     # Population mobility equilibrium condition;
     num=((a.^alpha).*(H.^(1-alpha)).*(dtradesh.^(-alpha./(sigma-1)))).^((sigma-1)./((sigma.*(1-alpha))-1));
     L_e=zeros(nobs,1);
     L_e = zeros(NN)
     L_e[1:div(NN, 2)] = num[1:div(NN, 2)] ./ sum(num[1:div(NN, 2)], dims=1)
     L_e[(div(NN, 2)+1):NN] = num[(div(NN, 2)+1):NN] ./ sum(num[(div(NN, 2)+1):NN], dims=1)
     L_e[1:div(NN, 2)] = L_e[1:div(NN, 2)] .* LLwest
     L_e[(div(NN, 2)+1):NN] = L_e[(div(NN, 2)+1):NN] .* LLeast


     # Convergence criterion;
     income_r=round.(income.*(10 .^6));
     expend_r=round.(expend.*(10 .^6));
     L_i_r=round.(L_i.*(10 .^6));
     L_e_r=round.(L_e.*(10 .^6));

     #disp([x max(abs(income_r-expend_r)) max(abs(L_e_r-L_i_r))]);

     # Update loop;
     if (income_r == expend_r) && (L_i_r == L_e_r)
        println(">>>> Convergence Achieved <<<<")
        x = 1000000
        converge = 1
      else
        w_e = w_i .* (expend ./ income).^(1 / (sigma - 1))
        w_i = (0.25 .* w_e) .+ (0.75 .* w_i)
        L_i = (0.25 .* L_e) .+ (0.75 .* L_i)
        # Normalization
        # Choose geometric mean wage in West as numeraire
        w_i[1:div(NN, 2)] = w_i[1:div(NN, 2)] ./ mean(w_i[1:div(NN, 2)])
        converge = 0
        x = x + 1
     end


    end

 # ****************************************************;
 # **** END LOOP TO SOLVE FOR POPULATION AND WAGES ****;
 # ****************************************************;

 xtic=time_ns();
 print(w_i)
 return w_i, L_i, tradesh, dtradesh, converge, xtic

end


