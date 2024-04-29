function RC_diffusion_factors_flux_control_elastic_1_continue(prefix, ...
        New_Terminating_time)
    % Continue running on previously-saved files
    % filename = strcat('../outputMats/', 'Elastic_Flag', num2str(Elastic_Flag), '_FluxTime_', num2str(flux), '_NewFH_', num2str(FHFlag), '_nuu_',  num2str(nuu), '_gamma_', num2str(gamma),...
    %                   '_pflag_', num2str(poreflag),'_c_', num2str(cc), '_factors_', ...
    %                   num2str(factors(1)), '_', num2str(factors(2)), '_',num2str(factors(3)),'.mat');
    
    filename = strcat('../outputMats/', prefix, '.mat');
    load(filename); 
    Terminating_time = New_Terminating_time; 
    NT = NT + it; 
    
    if ~exist('Elastic_Flag', 'var')
        Elastic_Flag = 0; 
    end
    
    % Main loop cannot be parallelized
    while it <= NT
        tryagaincount = 0;
        enter = 1;

        %% guess
        DX = (x(2)-x(1))^2;
        
        % Prescribed pc at the center of the injection
        if t + dt > 2144.2
            pcg_center = -191216.9;
        elseif t + dt > 1400
            pcg_center = -4544.5*(t+dt)+9553100;
        elseif t + dt > 1250
            pcg_center = 3190800;            
        else
            pcg_center = -2.073*(t + dt)^2+5144*(t+dt); % + 1.912e5; 
        end
        
        % Resolve the injection process, maximum increment in one step
        dtmax = 150 / 100;
        
        % Adapt injected mass accordingly
        % If Poroelastic
        if Elastic_Flag ~= 1 
            if it>2
                    % pcg = (- pave - (2*(bfs - bns)./(bfp + bnp).*siyy + 2*dphi./(phi*(bfp + bnp)) - 2*INjectmass(t+dt)*INprofile./(rhof0*phi*(bfp + bnp)) - KA*KD - KA*dt*(pave - pc) - KB*KDd - KB*dt*( pcpn -2*pc + pcmn + pavepn - 2*pave + pavemn)/DX ));
                    pcg_temp = (- pave - (2*(bfs - bns)./(bfp + bnp).*siyy + 2*dphi./(phi*(bfp + bnp)) - KA*KD - KA*dt*(pave - pc) - KB*KDd - KB*dt*( pcpn -2*pc + pcmn + pavepn - 2*pave + pavemn)/DX ));
                    % INjectma = (pcg_center - (pcg_temp(size(pcg_temp,1)/2) + pcg_temp(size(pcg_temp,1)/2+1))/2) ...
                    %     * (rhof0*phi*(bfp + bnp)) / max(INprofile);
                    pcg = pcg_temp + INjectmass(t+dt)*INprofile./(rhof0*phi*(bfp + bnp));
                
            else
                % INjectma = (pcg_center - (pc(size(pc,1)/2) + pc(size(pc,1)/2+1))/2) ...
                %     * (rhof0*phi*(bfp + bnp)) / max(INprofile);
                pcg = pc + INjectmass(t+dt)*INprofile./(rhof0*phi*(bfp + bnp));
            end
        
        % Elastic case
        else
            RHS = pcp - (bfs - bns)./(bfp + bnp) .* (siyy - siyyp) ...
                      - (dphi - dphip) ./ (phi*(bfp + bnp)) ...
                      + (INjectmass(t+dt) - INjectmass(t)) * INprofile./(rhof0*phi*(bfp + bnp)) ...
                      ;
            M = spdiags(1 + dt*(2*KB/DX)*ones(length(x),1),0,length(x),length(x)) ...
                + spdiags(-dt * KB/DX * ones(length(x),2),[-1 1],length(x),length(x));
            M(1,length(x)) = -dt*(KB)/DX;
            M(length(x),1) = -dt*(KB)/DX;
            pcg = M \ RHS; 
            % pavepn = [pave(end);pave(1:end-1)];
            % pavemn = [pave(2:end);pave(1)];
            % pavepp = [pavep(end);pavep(1:end-1)];
            % pavemp = [pavep(2:end);pavep(1)];
            % pcpp = [pcp(end);pcp(1:end-1)];
            % pcmp = [pcp(2:end);pcp(1)];
            % 
            % LK = - KB*KDdp - 0.5*KB*dt*(pcpp - 2*pcp + pcmp + pavepp - 2*pavep + pavemp + pavepn - 2*pave + pavemn)/DX;
            % 
            % M = spdiags(1+0.5*dt*(2*KB/DX)*ones(length(x),1),0,length(x),length(x)) + spdiags(-0.5*dt*KB/DX*ones(length(x),2),[-1 1],length(x),length(x));
            % M(1,length(x)) = -0.5*dt*(KB)/DX;
            % M(length(x),1) = -0.5*dt*(KB)/DX;
            % 
            % pcg = M\(- pave - (2*(bfs - bns)./(bfp + bnp).*siyy + 2*dphi./(phi*(bfp + bnp)) - INjectmass(t+dt)*INprofile./(rhof0*phi*(bfp + bnp)) + LK ));
        end
        
        dyg = dy + dt*(dy-dyp)/dtp;
        Vg =  V + dt*(V-Vp)/dtp;
        Vmg = 0.5*(V + Vg);
        dxg = dx + dt*Vmg;

        pavep = pave;
        sigrp = sigr;
        sigrnp = sigrn;
        pcp = pc;
        dyp = dy;
        dxp = dx;

        % Update siyyp and dphip
        siyyp = siyy; 
        dphip = dphi;

        thetap = theta;
        KDp = KD;
        KDdp = KDd;
        Vp = V;

        while tryagain || enter
            counter = counter + 1;

            % If this is the first attempt
            if enter == 1
                F =  fftshift(fft(dxg))/N;
                Fy = fftshift(fft(dyg))/N;
                Fp = fftshift(fft(pcg))/N;

            elseif tryagaincount > 1
                pcg = 0.5*(pc + pcp);
                dxg = 0.5*(dx + dxp);
                dyg = 0.5*(dy + dyp);
                F =  fftshift(fft(dxg))/N;
                Fy = fftshift(fft(dyg))/N;
                Fp = fftshift(fft(pcg))/N;
                Vg = 0.5*(Vp+V);
                Vmg = 0.5*(Vp + Vg);
            elseif tryagaincount == 1
                %first iteration is treated as a special case where time-step is not decreased
                % in attempt to optain a better guess
                pcg = pc;
                dxg = dx;
                dyg = dy;
                F =  fftshift(fft(dxg))/N;
                Fy = fftshift(fft(dyg))/N;
                Fp = fftshift(fft(pcg))/N;
                Vg = V;
                Vmg = 0.5*(Vp + Vg);

            end
            enter = 0;


            if AL == 1
                thetag = thetap.*exp(-Vmg*dt/L) + L./Vmg.*(1 - exp(-Vmg*dt/L));
            else
                thetag = L./Vmg.*(Vmg.*thetap/L).^(exp(-Vmg*dt/L));
            end
            %% Original way of doing convolution
            % Interpolate to get Fourier of dx at (t+dt+TT), namely DDT
            %DDTcor = (t + dt - tup)*(F-DD)./(-TT + t + dt - tup);
            %DDT = DD + DDTcor;

            %DDyTcor = (t + dt - tup)*(Fy-DDy)./(-TT + t + dt - tup);
            %DDyT = DDy + DDyTcor;

            %PPTcor = (t + dt - tup)*(Fp-PP)./(-TT + t + dt - tup);
            %PPT = PP + PPTcor;


            %taucon   = -TT(:,end).*K1(:,end).*DDT(:,end) +  0.5*sum((TT(:,2:end) - TT(:,1:end-1) ).*(K1(:,2:end).*DDT(:,2:end) + K1(:,1:end-1).*DDT(:,1:end-1)),2);
            %sigcon   = -TT(:,end).*K2(:,end).*DDT(:,end) +  0.5*sum((TT(:,2:end) - TT(:,1:end-1) ).*(K2(:,2:end).*DDT(:,2:end) + K2(:,1:end-1).*DDT(:,1:end-1)),2);
            %sigcony  = -TT(:,end).*K2(:,end).*DDyT(:,end) + 0.5*sum((TT(:,2:end) - TT(:,1:end-1) ).*(K2(:,2:end).*DDyT(:,2:end) + K2(:,1:end-1).*DDyT(:,1:end-1)),2);
            %sigconp  = -TT(:,end).*K2(:,end).*PPT(:,end) +  0.5*sum((TT(:,2:end) - TT(:,1:end-1) ).*(K2(:,2:end).*PPT(:,2:end) + K2(:,1:end-1).*PPT(:,1:end-1)),2);
            %sigyycon = -TT(:,end).*K1(:,end).*PPT(:,end) +  0.5*sum((TT(:,2:end) - TT(:,1:end-1) ).*(K1(:,2:end).*PPT(:,2:end) + K1(:,1:end-1).*PPT(:,1:end-1)),2);
            %sigyycony= -TT(:,end).*K1(:,end).*DDyT(:,end) + 0.5*sum((TT(:,2:end) - TT(:,1:end-1) ).*(K1(:,2:end).*DDyT(:,2:end) + K1(:,1:end-1).*DDyT(:,1:end-1)),2);
            
            %% Doing convolution is needed for poroelasticity, not for elasticity
            if Elastic_Flag ~= 1
                % Parallel way of doing convolution
                DDTcor = (t + dt - tup)*(F-DD)./(-TT + t + dt - tup);
                DD_var(:,:,1) = DD + DDTcor;
    
                DDyTcor = (t + dt - tup)*(Fy-DDy)./(-TT + t + dt - tup);
                DD_var(:,:,2) = DDy + DDyTcor;
    
                PPTcor = (t + dt - tup)*(Fp-PP)./(-TT + t + dt - tup);
                DD_var(:,:,3) = PP + PPTcor;
    
    
                for i = 1:1:6
                    con_var(:,i) = - TT_K(:,ind_0(i)) .* DD_var(:,end,ind_1(i))...
                        + 0.5*sum(TT_diff .* ...
                        (K_var(:,2:end, ind_2(i)).*DD_var(:,2:end,ind_3(i))...
                        + K_var(:,1:end-1,ind_4(i)).*DD_var(:,1:end-1,ind_5(i))), 2);
                end
                
                % If fully poroelastic
                if Elastic_Flag == 0
                    % Slightly modified to match variable names
                    taur =  real(ifft(ifftshift(( akdhat1.*(con_var(:,1)+F)))));
                    sigr =  real(ifft(ifftshift(( akdhat2.*(con_var(:,2)+F)  + akdhat3.*con_var(:,4) + akdhat5.*(con_var(:,3)+Fy))))); %!
                    sigrn = real(ifft(ifftshift((-akdhat2.*(con_var(:,2)+F)  + akdhat3.*con_var(:,4) + akdhat5.*(con_var(:,3)+Fy))))); %!
                    siyy =  -real(ifft(ifftshift(( akdhat4.*con_var(:,5) + akdhat1.*(con_var(:,6)+Fy)))));
                
                % If the bulk is only permeable, not poroelastic (Elastic_Flag == 2)
                else
                    taur =  real(ifft(ifftshift(( akdhat1 .* F))));
                    sigr =  real(ifft(ifftshift((akdhat3.*con_var(:,4)))));
                    sigrn = sigr; 
                    siyy =  -real(ifft(ifftshift((akdhat1 .* Fy))));
                end

            else % Elastic case
                % Slightly modified to match variable names
                taur =  real(ifft(ifftshift(( akdhat1 .* F))));
                sigr =  pcg; %!
                sigrn = pcg; %!
                siyy =  -real(ifft(ifftshift((akdhat1 .* Fy))));
            end

            tau = tau0 + taur;

            if poreflag == 1
                si  =  si0 - sigr + siyy;
            elseif poreflag == 2
                si  =  si0 - sigrn + siyy;
            elseif poreflag == 3
                si  =  si0 - 0.5*(pcg + 0.5*(sigr + sigrn)) + siyy;
            elseif poreflag == 4
                si  =  si0 - pcg + siyy;
            elseif poreflag == 5
                sig = si0 + siyy;
                si = (sig./(pcg./(2*sig) + sigrn./(4*sig) + sigr./(4*sig) + 1)*tol + 2*(pcg-sigrn).*(pcg-sigr) )./(tol + log( ((sig-sigr)./(sig-pcg))).*(pcg-sigrn) + log((sig-sigrn)./(sig-pcg)).*(pcg-sigr) ) ;
            elseif poreflag == 6                
                % Assumes slip happens at minimum effective normal stress
                si = si0 + siyy - max([pcg'; sigr'; sigrn'])';
            end

            Isineg = si <= 0;
            si(Isineg) = 1;
            tau(Isineg) = 0;

            % If Vg large enough, use implicit method to determine V
            if max(Vg)/Vr >= Vthres % && rem(i,1)==0
                I = (Vg/Vr < Vthres) & (d0 == 1);
                II = find(Vg/Vr >= Vthres);
                if FHFlag == 0
                    f = @(x)abs((x/Vr - exp( ( (tau(II)-eta*x)./(si(II)) - fr - b(II).*log(thetag(II)/(L/Vr)) )./a(II) ))) ;%- const*([0;Vp(1:end-2) - 2*Vp(2:end-1) + Vp(3:end);0])/dt^2;
                    % f = @(x) sum(abs((x/Vr - exp( ( (tau(II)-eta*x)./(si(II)) - fr - b(II).*log(thetag(II)/(L/Vr)) )./a(II) )))) ;%- const*([0;Vp(1:end-2) - 2*Vp(2:end-1) + Vp(3:end);0])/dt^2
                else
                    f = @(x)abs(x/Vr - 2 * sinh((((tau(II) - eta * x) ./ (si(II)) - fw) .* (1 + L ./ thetag(II) ./ Vw) + fw) ./ a(II)) ./ exp((fr + b(II) .* log(Vr * thetag(II) / L)) ./ a(II)));
                    % f = @(x) sum(abs(x/Vr - 2 * sinh((((tau(II) - eta * x) ./ (si(II)) - fw) .* (1 + L ./ thetag(II) ./ Vw) + fw) ./ a(II)) ./ exp((fr + b(II) .* log(Vr * thetag(II) / L)) ./ a(II))));
                end
                VT = zeros(15,length(II));
                VpT = Vg(II);
                fact = 1 + 0.03*((linspace(-1,1,15)).^5) ;%linspace(0.98,1.03,length(VT(:,1)));

                % Can be parrelized
                for iII = 1:length(VT(:,1))
                    VT(iII,:) = f(fact(iII)*VpT);
                end
                if length(VT(1,:))>1
                    [~,imin] = min(VT,[],1);
                else
                    [~,imin] = min(VT,[],1);
                end
                V(II) = VpT.*fact(imin)';
                V(I) =   (Vr*exp( ( (tau(I)-eta*Vg(I))./(si(I)) - fr - b(I).*log(thetag(I)/(L/Vr)))./a(I)) );
                % options = optimset('TolFun', 1.e-14, 'TolX', 1.e-14);
                % [VpT, FpT000] = fminsearch(f, Vg(II), options); 
                % V(II) = VpT; 

                if FHFlag == 0
                    % Regularized friction law
                    V(I) = 2 * Vr * sinh((tau(I) - eta * Vg(I)) ./ (si(I)) ./ a(I)) ./ exp((fr + b(I) .* log(Vr * thetag(I) / L)) ./ a(I));
                else
                    % Regularized friction law with flash heating
                    V(I) = 2 * Vr * sinh((((tau(I) - eta * Vg(I)) ./ (si(I)) - fw) .* (1 + L ./ thetag(I) ./ Vw) + fw) ./ a(I)) ./ exp((fr + b(I) .* log(Vr * thetag(I) / L)) ./ a(I));
                end
            % If Vg small, use explicit scheme
            else
                %V = Vr*exp( ( (tau-eta*Vg)./(si) - fr - b.*log(thetag/(L/Vr)))./a) ;
                if FHFlag == 0
                    % Regularized friction law
                    V = 2 * Vr * sinh((tau - eta * Vg) ./ (si) ./ a) ./ exp((fr + b .* log(Vr * thetag / L)) ./ a);
                else
                    % Regularized friction law with flash heating
                    V = 2 * Vr * sinh((((tau - eta * Vg) ./ (si) - fw ) .* (1 + L ./ thetag ./ Vw) + fw) ./ a) ./ exp((fr + b .* log(Vr * thetag / L)) ./ a);
                end
            end
            V = V.*d0 + Vo.*d00;
            Vm = 0.5*(Vp + V);

            % Aging/Slipping law for state variable update
            if AL == 1
                theta = thetap.*exp(-Vm*dt/L) + L./Vm.*(1 - exp(-Vm*dt/L));
            else
                theta = L./Vm.*(Vm.*thetap/L).^(exp(-Vm*dt/L));
            end
            dphi = dphi0 - Gamma*log(Vr.*theta/L);
            pave = 0.5*(sigr + sigrn);
            pavepn = [pave(end);pave(1:end-1)];
            pavemn = [pave(2:end);pave(1)];
            pavepp = [pavep(end);pavep(1:end-1)];
            pavemp = [pavep(2:end);pavep(1)];
            pcpp = [pcp(end);pcp(1:end-1)];
            pcmp = [pcp(2:end);pcp(1)];

            LK = - KB*KDdp - 0.5*KB*dt*(pcpp - 2*pcp + pcmp + pavepp - 2*pavep + pavemp + pavepn - 2*pave + pavemn)/DX;

            % Poroelastic case
            if Elastic_Flag ~= 1
                M = spdiags(1+0.5*dt*(KA+2*KB/DX)*ones(length(x),1),0,length(x),length(x)) + spdiags(-0.5*dt*KB/DX*ones(length(x),2),[-1 1],length(x),length(x));
                M(1,length(x)) = -0.5*dt*(KB)/DX;
                M(length(x),1) = -0.5*dt*(KB)/DX;
            % Elastic case
            else
                M = spdiags(1+dt*(2*KB/DX)*ones(length(x),1),0,length(x),length(x)) ...
                    + spdiags(-dt*KB/DX*ones(length(x),2),[-1 1],length(x),length(x));
            end


            
            % Elastic case
            if Elastic_Flag ~= 1
                pc = M\(- pave - (2*(bfs - bns)./(bfp + bnp).*siyy + 2*dphi./(phi*(bfp + bnp)) - INjectmass(t+dt)*INprofile./(rhof0*phi*(bfp + bnp)) - KA*KDp - 0.5*KA*dt*(pave + pavep - pcp) + LK ));
            else
                RHS = pcp - (bfs - bns)./(bfp + bnp) .* (siyy - siyyp) ...
                      - (dphi - dphip) ./ (phi*(bfp + bnp)) ...
                      + (INjectmass(t+dt) - INjectmass(t)) * INprofile./(rhof0*phi*(bfp + bnp)) ...
                      ;
                pc = M \ RHS; 
            end

            pm = 0.5*(pc + 0.5*(sigr + sigrn));
            dy = 2*epsi*(phi/(1-phi)*bnp - bgp).*(pm - (phi/(1-phi)*bns - bgs)./(phi/(1-phi)*bnp - bgp).*siyy) + 2*epsi*dphi/(1-phi);
            dx = dxp + dt*Vm;

            accuracyrel = norm(pc-pcg,1)/(norm(pc,1));
            [accuracyabs,~] = max(abs((pc-pcg))./(a.*si0));

            pcpn = [pc(end);pc(1:end-1)];
            pcmn = [pc(2:end);pc(1)];
            
            % IF poroelastic, no cross-fault fluid motion is allowed
            if Elastic_Flag ~= 1
                KD = KDp + 0.5*dt*(pave-pc + pavep-pcp);
            
            % Elastic
            else
                KD(:) = 0.;
                sigr = pc;
                sigrn = pc; 
                pave = 0.5*(sigr + sigrn);
                pavepn = [pave(end);pave(1:end-1)];
                pavemn = [pave(2:end);pave(1)];
                pavepp = [pavep(end);pavep(1:end-1)];
                pavemp = [pavep(2:end);pavep(1)];
                pcpp = [pcp(end);pcp(1:end-1)];
                pcmp = [pcp(2:end);pcp(1)];
            end 
            
            KDd = KDdp + 0.5*dt*( pcpn -2*pc + pcmn + pcpp - 2*pcp + pcmp + pavepp - 2*pavep + pavemp + pavepn - 2*pave + pavemn)/DX;
            
            if (tryagaincount < miniter || max(abs((pc-pcg))./(a.*si0)) > 0.1*frac || norm(pc-pcg,1)/(norm(pc,1)+1) > 0.1*frac) && tryagaincount < 10

                tryagain = 1;
                tryagaincount = tryagaincount + 1;
                if tryagaincount==1
                    % First iteration refines the prediction and does not reduce the time-step
                    % dt = dt;
                else
                    % If accuracy metrics are violated then time-step is refined
                    dt = dt/2;
                end

            else
                tryagain = 0;
                tryagaincount = 0;
                cumutolvio = cumutolvio + dt*accuracyabs;
                if accuracyabs > maxtolvio
                    maxtolvio = accuracyabs;
                    maxtolviotime = t;
                end
            end


            %% Update kernels
            % Update happens after this time-step has been fully resolved
            if tryagain == 0

                if rem(it,NSplot)==0 || it == 2 || (runnerplot > 1 && max(V./Vsave(:,runnerplot - 1)) > 10)
                    Vsave(:,runnerplot) = V;
                    dsave(:,runnerplot) = dx;
                    psave(:,runnerplot) = pm;
                    tsaveplot(runnerplot) = t;
                    dysave(:,runnerplot) = dy;
                    sigrsave(:,runnerplot) = sigr;
                    pcsave(:,runnerplot) = pc;
                    sisave(:,runnerplot) = si;
                    tauS(:,runnerplot) = tau;
                    thetasave(:,runnerplot) = theta;
                    InjectMaSave(runnerplot) = INjectmass(t+dt); 
                    if plotflag == 1
                        f1 = figure(1);
                        sgtitle(strcat('\epsilon = ',' ',num2str(epsi),' --- ','time =',' ',num2str(t/(24*60*60)),' ',' days'))
                        subplot(3,2,1)
                        semilogy(x,V,'k-')
                        xlim([min(x) max(x)])
                        %hold on
                        title('V (slip speed) [m/s]')
                        subplot(3,2,5)

                        plot(x,dy*1.0e3,'r-')
                        xlim([min(x) max(x)])
                        ylim([-3 3]*1.0e-2);

                        title('d_y (normal displacements) [mm]')
                        subplot(3,2,4)
                        plot(x,tau/1.0e6,'r-')
                        xlim([min(x) max(x)])
                        ylim([tau0/1.0e6-3 tau0/1.0e6+3])
                        %hold on
                        title('tau (total shear stress) [MPa]')

                        subplot(3,2,3)
                        plot(x,dx,'b-')

                        xlim([min(x) max(x)])

                        title('d_x (cumulative slip) [m]')
                        subplot(3,2,2)

                        plot(x,pm/1.0e6,'-')
                        hold on
                        plot(x,sigr/1.0e6,'-')
                        plot(x,sigrn/1.0e6,'-')
                        plot(x,pc/1.0e6,'-')
                        xlim([min(x) max(x)])
                        ylim([-1 2.5])
                        hold off
                        legend('p_m','p^+','p^-','p_c');

                        %    hold on
                        title('pore-pressure [MPa]')
                        %ylim([-0.05 0.05])
                        subplot(3,2,6)

                        plot(x,si/1.0e6,'b-')
                        xlim([min(x) max(x)])
                        %hold on
                        title('sig (total effective normal stress) [MPa]')
                        ylim([si0/1.0e6-0.2 si0/1.0e6+0.2])
                        set(gcf,'color','w');
                        drawnow
                        if gifflag
                            if firstframe
                                gif(GIFNAME,'DelayTime',0.25,'frame',f1,'nodither')
                                firstframe = 0;
                            else
                                gif
                            end
                        end
                    end
                    runnerplot = runnerplot + 1;
                    %it
                end


                diftime = min([(1/(c*(2*pi)^2/(x(2)-x(1))^2)) 1/(c)*((kappa*epsi)/(kappac))^2]);
                NS = floor(min([1+CNS*((diftime/dt)) 1+0.05*min((a*si0./abs(pc-pcNSp))) MAXNS]));
                if t+dt-tup < alpmin*diftime
                    NS = MAXNS;
                end

                % If iteration steps exceed max number of steps before kernel
                % update
                if Kit + NS < it + 1 || it == 2
                    Kit = it;
                    pcNSp = pc;
                    tsave(runner) = t+dt;

                    % Save for kernel updates:
                    FR(:,runner) = real(F);
                    FI(:,runner) = imag(F);
                    FRy(:,runner) = real(Fy);
                    FIy(:,runner) = imag(Fy);
                    FRp(:,runner) = real(Fp);
                    FIp(:,runner) = imag(Fp);

                    if runner > 1

                        % Calculate TTI outside the loop to save time
                        TTI = TT' + t + dt;

                        % This loop can be parallelized
                        for ixx = 1:length(x)
                            ITT = find(tsave(1:runner) < (TT(ixx,1)+t),1,'last' )  ;
                            if isempty(ITT)==1
                                ITT = 1;
                            end

                            samplepoints = ITT:runner;

                            % Calculate TTI outside the loop can save time
                            %TTI = TT(ixx,:)'+ t + dt;

                            FRI = FR(ixx,samplepoints)';
                            FII = FI(ixx,samplepoints)';

                            % Change TTI to TTI(:,ixx)
                            DD(ixx,:) = interp1q(tsave(samplepoints)',FRI,TTI(:,ixx)) ...
                                + 1i*interp1q(tsave(samplepoints)',FII,TTI(:,ixx)) ;

                            % Need to move the following line outside the loop
                            % for parallelization, same for DDy and PP
                            %DD(ixx,isnan(DD(ixx,:))) = 0;

                            %VV(ixx,:) = [diff(DD(ixx,:))./diff(TT(ixx,:)) 0];
                            %VV(ixx,isnan(VV(ixx,:))) = 0;

                            % Not important if gamma = 0
                            if Gamma ~= 0
                                FRI = FRy(ixx,samplepoints)';
                                FII = FIy(ixx,samplepoints)';
                                DDy(ixx,:) = interp1q(tsave(samplepoints)',FRI,TTI(:,ixx))...
                                    + 1i*interp1q(tsave(samplepoints)',FII,TTI(:,ixx)) ;
                                %DDy(ixx,isnan(DDy(ixx,:))) = 0;
                                %VVy(ixx,:) = [diff(DDy(ixx,:))./diff(TT(ixx,:)) 0];
                                %VVy(ixx,isnan(VVy(ixx,:))) = 0;
                            end

                            FRI = FRp(ixx,samplepoints)';
                            FII = FIp(ixx,samplepoints)';
                            PP(ixx,:) = interp1q(tsave(samplepoints)',FRI,TTI(:,ixx))...
                                + 1i*interp1q(tsave(samplepoints)',FII,TTI(:,ixx)) ;
                            %PP(ixx,isnan(PP(ixx,:))) = 0;
                            %PV(ixx,:) = [diff(PP(ixx,:))./diff(TT(ixx,:)) 0];
                            %PV(ixx,isnan(PV(ixx,:))) = 0;

                            % Keep record of last kernel update time (Moved out of this loop???)
                            %tup = t + dt;
                        end

                        % Eliminate NaNs
                        DD(isnan(DD)) = 0;
                        DDy(isnan(DDy)) = 0;
                        PP(isnan(PP)) = 0;

                        % Keep record of last kernel update time
                        tup = t + dt;
                    end

                    % Runner keeps record of kernel update steps
                    runner = runner + 1;

                end

                if isnan(sum(V))
                    % Truncate the zero endings
                    Vsave = Vsave(:,1:runnerplot - 1);
                    dsave = dsave(:,1:runnerplot - 1);
                    psave = psave(:,1:runnerplot - 1);
                    tsaveplot = tsaveplot(:,1:runnerplot - 1);
                    dysave = dysave(:,1:runnerplot - 1);
                    sigrsave = sigrsave(:,1:runnerplot - 1);
                    pcsave = pcsave(:,1:runnerplot - 1);
                    sisave = sisave(:,1:runnerplot - 1);
                    tauS = tauS(:,1:runnerplot - 1);
                    thetasave = thetasave(:, 1:runnerplot - 1);
                    InjectMaSave = InjectMaSave(:, 1:runnerplot - 1); 

                    % Filename reflects fract number and parallelization
                    filename = strcat('../outputMats/', prefix, '_continue_', num2str(New_Terminating_time), '.mat');
                    % Record excuting time of the program
                    %t1 = cputime - t0;
                    t1 = toc(tstart);
                    disp('Finished!');
                    save(filename);

                    % Write changable parameters into a '.txt' file
                    txtname = strcat('../outputMats/', prefix, '_continue_', num2str(New_Terminating_time), '.txt');

                    fileID = fopen(txtname, 'w');

                    fprintf(fileID, '%25s', 'Filename: '); fprintf(fileID, filename);
                    fprintf(fileID, '\n%25s', 'NT: '); fprintf(fileID, num2str(NT));
                    fprintf(fileID, '\n%25s', 'lhs, rhs: '); 
                    fprintf(fileID, strcat(num2str(lhs), ', ', num2str(rhs))); 
                    fprintf(fileID, '\n%25s', 'Gamma: '); fprintf(fileID, num2str(Gamma));
                    fprintf(fileID, '\n%25s', 'a, b: '); 
                    fprintf(fileID, strcat(num2str(a(1)), ', ', num2str(b(1)))); 
                    fprintf(fileID, '\n%25s', 'Terminated Slip rate: '); fprintf(fileID, num2str(Terminating_slip_rate));
                    fprintf(fileID, '\n%25s', 'Total Wallclock Time (min): '); fprintf(fileID, num2str(t1/60));
                    fprintf(fileID, '\n%25s', 'Number of iterations:'); fprintf(fileID, num2str(it));
                    fprintf(fileID, '\n%25s', 'kappac:'); fprintf(fileID, num2str(kappac));
                    fprintf(fileID, '\n%25s', 'kappacx:'); fprintf(fileID, num2str(kappacx));
                    fprintf(fileID, '\n%25s', 'c:'); fprintf(fileID, num2str(c));
                    fclose(fileID);
                    disp('NaN detected in V!');
                    disp(strcat('Gamma = ', num2str(Gamma)));
                    return;
                end
                t = t + dt;
                dtp = dt;

                % Reset dt for the next time step
                dt = min([frac*min([L./max(V)]), dtmax]);
                if dt > maxdtfac*dtp
                    dt = maxdtfac*dtp;
                end
            end
        end
        
        % Stop after 12 days if no injection happens
        if t > Terminating_time
            break;
        end
        
        % Stop after 0.1 m/s
        %if runnerplot > 2 &&  Vsave(64,runnerplot-2) > Terminating_slip_rate && V(64) <= Terminating_slip_rate
        %    break;
        %end
        it = it + 1;
    end

    % Truncate the zero endings
    Vsave = Vsave(:,1:runnerplot - 1);
    dsave = dsave(:,1:runnerplot - 1);
    psave = psave(:,1:runnerplot - 1);
    tsaveplot = tsaveplot(:,1:runnerplot - 1);
    dysave = dysave(:,1:runnerplot - 1);
    sigrsave = sigrsave(:,1:runnerplot - 1);
    pcsave = pcsave(:,1:runnerplot - 1);
    sisave = sisave(:,1:runnerplot - 1);
    tauS = tauS(:,1:runnerplot - 1);
    thetasave = thetasave(:, 1:runnerplot - 1);
    InjectMaSave = InjectMaSave(:, 1:runnerplot - 1);
    
    % Filename reflects fract number and parallelization
    % Filename reflects fract number and parallelization
    filename = strcat('../outputMats/', prefix, '_continue_', num2str(New_Terminating_time), '.mat');
    
    % Record excuting time of the program
    % t1 = cputime - t0;
    t1 = toc(tstart);
    disp('Finished');
    save(filename);
    
    % Write changable parameters into a '.txt' file
    txtname = strcat('../outputMats/', prefix, '_continue_', num2str(New_Terminating_time), '.txt');

    fileID = fopen(txtname, 'w');
    
    fprintf(fileID, '%25s', 'Filename: '); fprintf(fileID, filename);
    fprintf(fileID, '\n%25s', 'NT: '); fprintf(fileID, num2str(NT));
    fprintf(fileID, '\n%25s', 'lhs, rhs: '); 
    fprintf(fileID, strcat(num2str(lhs), ', ', num2str(rhs))); 
    fprintf(fileID, '\n%25s', 'Gamma: '); fprintf(fileID, num2str(Gamma));
    fprintf(fileID, '\n%25s', 'a, b: '); 
    fprintf(fileID, strcat(num2str(a(1)), ', ', num2str(b(1)))); 
    fprintf(fileID, '\n%25s', 'Terminated Slip rate: '); fprintf(fileID, num2str(Terminating_slip_rate));
    fprintf(fileID, '\n%25s', 'Total Wallclock Time (min): '); fprintf(fileID, num2str(t1/60));
    fprintf(fileID, '\n%25s', 'Number of iterations:'); fprintf(fileID, num2str(it));
    fprintf(fileID, '\n%25s', 'kappac:'); fprintf(fileID, num2str(kappac));
    fprintf(fileID, '\n%25s', 'kappacx:'); fprintf(fileID, num2str(kappacx));
    fprintf(fileID, '\n%25s', 'c:'); fprintf(fileID, num2str(c));
    fclose(fileID);
    
    disp(strcat("Time cost: ", num2str(t1)));
end

