function Rc_diffusion_factors_flux_control_elastic_1_0826(Nuu, Nu, ... 
                    Gamma, cc, ...
                    FHFlag, poreflag, factors, flux, Elastic_Flag, terminating_time, ...
                    initial_V_ratio)
    %% factors: 
    % diffusivity factor for [kappacx, kappacy and c] of the bulk 
    %% Problem setup 
    % IMPORTANT: This code modified to replicate Stacy's paper results with
    % INTERMEDIATE friction parameters
    % this code set up an injection problem into a fault into a region with a width of 30m. Fluid is injected at
    % high rate for 1 h. This causes a delayed nucleation around 3.2 days if
    % time-integration is suffuciently accurate. The shear zone has across
    % fault permeability/mobilit that is 7 orders of magnitude smaller than the
    % along fault permeability.
    % poreflag -- 3 average, 6 maximum
    %% Parallelization setup
    % Several "parfor" loops are used, 
    % the number of workers for parfor is set to be default value (<= 12 on my computer)
    % it can also be explicitly set via
    % "parfor (i = 1:1:n, M)"
    % where M, a non-negative integer specifies maximum number of workers

    % Start clock
    % t0 = cputime;
    tstart = tic;
    
    %% IMPORTANT: ELASTIC FLAG: 1, elastic; 0, normal; 
    % Elastic by default uses nu_u
    % Elastic_Flag = 0; 
    
    % Terminating slip rate and simulating time
    Terminating_slip_rate = 1.0e-1;
    baseFlux = 1.0e-4; 
    Terminating_time = terminating_time; % baseFlux * 2020 / flux;
    if Terminating_time == 0
        Terminating_time = baseFlux * 2020 / flux * 4.;
    end
    
    % Terminating time if in_mass = 0, 12 days
    % if in_mass == 0
    %     Terminating_time = 10000;
    % end
    %% IMPORTNANT: TIME-STEPPING SCALING: frac
    fract = 4; % used in frac to scale time-steps You may want to pick a smaller value
    % when trying things out to speed things up
    miniter = 1; % minimum number of iterations. For the updates scheme it is best to make this 1.

    % frac: factor scaling the adaptive time-step
    % you generally want to try running simulations where this is smaller to
    % check if everything is converging. How small this needs to be depends on
    % the problem you are running.

    %if NS, which gets updated and printed out every time-step is 1 that is an indication that frac should be smaller
    frac = 0.5/2^(fract); %1/2^3 is at the upper limit to what is acceptable, larger frac may cause instability
    % note that frac also scales the error tolerance for pore-pressure

    %%
    rng(1); %set random seed

    %% FIGURE
    % the script plots some fault fields as it runs. This can be deactivated by
    % setting plotflag = 0

    plotflag = 0;

    if plotflag
        f1 = figure(1);
        set(gcf, 'Position', [0 0 1200 1500]); %this may need to be changed if the figure doesn't show up
    end
    %% MAKE A GIF
    %if you want to ex make a gif set gifflag = 1
    gifflag = 0;
    firstframe = 1;
    GIFNAME = 'TEST_injection.gif';

    %% SELECT A PORE-PRESSURE MODEL
    %poreflag = 1; % means when you compute the effective normal stress you use
    %the porepressure on one side of the shear zone (y positive side). SAME AS
    %JMPS!

    %poreflag = 2; %uses the y negative side

    %poreflag = 3; %uses the average porepressure (this will be used by Heimisson, Lapusta and Rudnicki)

    %poreflag = 4; %uses the central pore - pressure pc

    %poreflag = 5; %true average pore pressure
    
    %poreflag = 6; %uses the maximum of p+, p-, pc

    %% SELECT FRICTION LAW
    %AL = 1 % select ageing law
    AL = 1; % or anything else will use slip-law (used in JMPS)
    %%
    % All units are standard SI

    % x cooridinates
    % x = linspace(-250, 250, 2^10);
    x = linspace(-312.5, 312.5, 1280);
    % x = linespace(-500, 500, 2^11);
    N = length(x);
    % KL: number of values used in computing convolution kernels
    % large KL is more accurate but slower, I've found values in the range between 2^6-2^8 to be sufficient
    KL = 2^10;
    %convolution goes to alp*[max diffusion time-scale]:
    alp = 20;
    %convolution starts at alpmin*[min diffusion time-scale]:
    alpmin =  1.0e-6;
    
    % Rate and State parameters
    b = 0.0160;
    a = 0.01125;
    % Vharacteristic state evolution distance
    L = 16.75e-6;
    % Initial state variable
    theta_0 = 2.38e12;
    
    %initial and reference stresses
    fr = 0.55;
    f0 = 0.5375;
    si0 = 4e6;
    tau0 = f0*si0;
    
    %Flash heating constants
    fw = 0.3;
    Vw = 1e-1;
    
    % Calculate initial slip rate
    Vr = 1.0e-6; % Reference slip rate
    Vo = 2 * Vr * sinh(tau0 / a / si0) / exp((fr + b * log(Vr * theta_0 / L)) / a); % Initial slip rate
    Vo = Vo * initial_V_ratio; 
    theta_0 = L / Vr * exp((a * log(2 * Vr / Vo * sinh(tau0 / a / si0)) - fr) / b);


    NT = 30000000;% max number of time-steps
    Vthres = 1.0e4; % threshold max(V)/Vo ratio at which an implicit step is taken
    % Vthres = 0.;

    %Vthres > 1.0e8, can cause spurious oscillations

    MAXNS = inf; %maximum number of timesteps needed for a kernel update
    % This limit may be reached if slip speed is high
    NS = 10; %gets changed, just to initialize and allowcate imperfectly
    NSplot = 10; % 10;
    %should here represent the average NS in the simulation.
    CNS = 1; % if larger than 1 the simulations may not fully resolve the shortest diffusion times, which typically depends on the
    %grid size, but will make them more efficient

    maxdtfac = 2^(miniter-1) + 0.1;
    if miniter == 0
        maxdtfac = 1.1;
    end
    Kit = 0;

    %defines a loading region, for compatison with JMPS we make this logical
    %mask larger then the domain
    d0 = (heaviside(-x + 200).*heaviside(x + 200))';
    % d0 = (heaviside(-x + 400).*heaviside(x + 400))';
    d00 = d0 == 0;


    %allocate slip
    dx = 0*d0;
    %allocate opening (only matters when dilatancy is included)
    dy = dx;
    %allocate slip zone central pressure
    pc = dx;
    %initialize slip speed
    V = dx + Vo;
    %innitialize Fourier coefficienents
    F = fftshift(fft(dx));
    Fs = 1./(x(2)-x(1));
    vn = Fs*(-N/2:1:N/2-1)/N;

    % Spatial frequency
    fre = vn;
    F = F/N;

    %material parameters
    nu = Nu;
    nuu = Nuu;
    B = 0.85;
    G = 10e9;
    rhof0 = 1.0e3; %reference fluid density kg/m^3

    %hydraulic diffusivity, by factor
    c = cc * factors(3);

    %kappac = 1.0e-17; %make the shear zone very impermeable to replicate JMPS
    %kappacx = 1.0e9*kappac;

    
    %(which assumes either fully impermeable or fully permeable)
    %computing mobility of host rock (dependent on c)
    alpB = 3*(nuu-nu)/(B*(1+nuu)*(1-2*nu));
    kappa = c/(2*G*(1+nu)*B / (3*alpB*(1-alpB*B)*(1-2*nu) ) );
    c = c*((1-nu)*(1-2*nuu) )/( (1-nuu)*(1-2*nu));
    
    % For elastic limit
    if Elastic_Flag == 1
        nuu = nu;
        c = 0.; 
    
    % Elastic bulk but with permeability
    elseif Elastic_Flag == 2
        nuu = nu; 
        alpB = 0.; 
        B = 0.; 
        kappa = 2.168791096126737e-19; % Computed with nu = 0.24, nuu = 0.35.
        c = cc * factors(3); 
    end
    
    %shear zone half thickness, doesn't really matter in comparing to JMPS
    %where the fault in completely impermeable
    epsi = 0.001;
    %shear wave velocity for radiaion damping
    cs = 3400;

    %radiation damping
    eta = G/(2*cs);

    %special compressibilities that only play a role when dilaticy is included
    bfp = 0.44e-9;
    bnp = 6.0e-9;
    bgp = 1/(50.0e9);

    bfs = 5*bfp/9;
    bns = 5*bnp/9;
    bgs = 5*bgp/9;
    %shear zone porosity
    phi = 0.068;

    % Make effective diffusivity 0.2 
    % alpha_x = kappacx / (phi * (bnp + bfp))
    
    % To recover 0.2, kappacx = 8.7584e-11
    kappacx = 8.7584e-11;
    
    L_nu = G * L / (b - a) / si0;
    % kappac = kappacx * factor * epsi^2 / L_nu^2;
    kappac = kappacx * 1.0e-9;

    % Scale kappacx and kappac
    kappacx = kappacx * factors(1);
    kappac = kappac * factors(2);
    
    %mildly rate-strengthening, + dx is just to generate a vector, this can be
    %spatially heterogeneous
    b = b + dx;
    a = a + dx;
    KA = 2*kappac./(epsi^2*phi.*(bfp + bnp));
    KB = kappacx./(phi.*(bfp + bnp));
    %% DILATANCY
    % Dilaticany coefficient Segall and Rice call this epsilon
    % JMPS code does not do dilaticany and this we put this to 0.
    gamma = Gamma;
    % gamma = 0.;
    %
    %% Injector
    %Gaussian
    %sigma = 50;
    %center = 0;
    %INprofile =  (1/(2*epsi*sigma*sqrt(2*pi))*exp(-0.5*((x-center)/sigma).^2))';
    %boxcar injector
    lhs = - 0.5;
    rhs = 0.5;
    %source term relative volume per unit length, should normalize to
    %1/2*epsi


    INprofile =  (heaviside(x - lhs).*heaviside(- x + rhs)/(2*epsi*(rhs-lhs)))';


    %%% NOTICE used to be volume and is now mass:
    % Injection over 1400s
    %in_mass = 3.386e6 * rhof0 * phi * (bfp + bnp) * 2 * epsi * (rhs - lhs);
    %in_rate = in_mass / 1400;
    %INjectmass = @(t) (in_rate) * t - ...
    %    (in_rate)*(t-1400)*heaviside(t-1400); % multiplies INprofile represents cumulative mass;
    %This is mathematically cleaner in the fluid mass balance
    %the source term is then INjectmass(t)*INprofile
    %load("./InjectMass.mat", "InjectMaSavesource", "tsaveplotsource");
    % INjectmass = @(t) interp1(tsaveplotsource, InjectMaSavesource, t, "spline");
    INjectmass = @(t) flux * t;
    %%
    % Initialize the state variable with small noise around steady state
    % NOTE that this is not the same noise as in the JMPS paper so you will see
    % qualitatively the same behavior but not exactly the same.
    % theta = L/Vo;
    
    theta = theta_0.*d0 + theta_0.*d00;
    dphi0 = gamma * log(Vr .* theta/L);
    % initialize time
    t = 0;
    % initialize time since last kernel update
    tup = 0;
    % make the first timestep something very small
    dt = 1.0e-5;
    % Previous time step value
    dtp = dt;
    % Initialize and integral kernel (only matters )
    KD = 0;
    KDd = 0;
    runner = 1;
    runnerplot = 1;
    tol = 1.0e-5;

    % Initialize a bunch of other things
    % Suffix "g" stands for guess

    % 
    % Pore pressure on positive side
    sigr = dx;

    % Pore pressure on negative side
    sigrn = dx;

    % Pore pressure on negative side, previous timestep
    sigrnp = dx;

    % Pore pressure on positive side, previous timestep
    sigrp = dx;

    % Guessed "change in total normal stress sig_yy"
    sigyyg = dx;

    % Guessed change in porosity IN THE SHEAR LAYER
    phig = dx;

    % Delta phi in the time step
    dphi = dx;

    % Change in total normal stress (sig_yy)
    siyy = dx;
    siyyp = dx; 

    % Opening, and previous step opening
    dy = dx;
    dyp = dy;

    % Averaged pore pressure with linear assumption, and previous step one
    pmp = dx;
    pm = dx;

    % Previous step center pressure
    pcp = dx;

    % Center pressure at the last kernel update
    pcNSp = dx;

    % Previous step slip
    dxp = dx;

    % Previous step state variable theta
    thetap = theta;

    % Pressure calculation variables ???
    KDp = KD;
    KDdp= KD;

    % Previous step V and dphi
    Vp = V;
    dphip = dx;

    % Average pressure, calculated by (sigr + sigrn)/2
    pave = dx;
    pavep = dx; 


    % Time history of a bunch of variables, at each kernel update

    % Time history of tauS
    tauS = zeros(length(x),NT/1000);

    % Time history of shear stress tau
    dsave = tauS;

    % Time history of slip rate V
    Vsave = tauS;

    % Time history of averaged pressure pm
    psave = tauS;

    % A few other saved variables for plot ---Shengduo
    dysave = tauS;
    sigrsave = tauS;
    sisave = tauS;
    pcsave = tauS;
    thetasave = tauS;

    % Real part of Fourier transform of slip, at each kernel update
    FR = tauS;

    % Imaginary part of Fourier transform of slip, at each kernel update
    FI = tauS;

    % Real part of Fourier transform of center pressure pc
    FRp = FR;

    % Imaginary part of Fourier transform of center pressure pc
    FIp = FI;

    % Real part of Fourier transform of opening dy, at each kernel update
    FRy = FR;

    % Imaginary part of Fourier transform of opening dy, at each kernel update
    FIy = FI;

    % Time of the saved data above
    tsave = zeros(1,NT/NSplot + 1);

    % Time of plot
    tsaveplot = zeros(1,NT/NSplot + 1);

    % Injection masses
    InjectMaSave = zeros(1, NT/NSplot + 1);


    %% Compute kernels
    FF = @(k)(kappac./(kappa*epsi*abs(k)));
    G1 = @(t,k)(- 2*(nuu-nu)/(1-nu)*c*k.^2.*(1 + FF(k)).*(1 + 1./(FF(k)-1) .* (FF(k).* (exp((FF(k).^2 - 1).*c*k.^2*t).*erfc(FF(k).*sqrt(c*k.^2*t)) - 1) + erf(sqrt(c*k.^2*t))) ));
    G2 = @(t,k)(- c*k.^2.* (1 + FF(k)).*( exp(-c*k.^2*t)./sqrt(pi*c*k.^2*t) - FF(k).*exp((FF(k).^2 - 1).*c*k.^2*t).*erfc(FF(k).*sqrt(c*k.^2*t)) )       );

    % Some variables for calculating time convolution

    % Integration nodes of numerical convolution
    TT = zeros(length(x),KL);
    TTp = TT;

    % Fourier transform of slip dx at (tup + TT)
    DD = TT;

    % Fourier transform of slip rate V, at (tup + TT)
    VV = TT;

    % Fourier transform of opening dy, at (tup + TT)
    DDy = TT;

    % Fourier transform of opening rate Vy, at (tup + TT)
    VVy = TT;

    % Fourier transform of center pressure pc, at (tup + TT)
    PP = TT;

    % Fourier transform of center pressure pc, at (tup + TT)
    PV = TT;

    % Kernels evaluated at (-TT)
    K1 = TT;
    K2 = TT;


    %% distretize the time-steps of the kernels (messier then it needs to be)
    % Can be parallelized since there's no interaction within the loop
    for i = 1:length(x)
        kv = fre(i)*(2*pi);

        if kv ~= 0
            ttt1 = logspace(log10(alpmin*min([1/(c*kv.^2) ((kappa*epsi)/(kappac))^2/c])),log10(alp*max([1/(c*kv.^2) ((kappa*epsi)/(kappac))^2/c])),KL/2);
            ttt2 = linspace((alpmin*min([1/(c*kv.^2) ((kappa*epsi)/(kappac))^2/c])),(alp*max([1/(c*kv.^2) ((kappa*epsi)/(kappac))^2/c])),KL/2+2);
            ttt2 = ttt2(2:end-1);
            ttt = sort([ttt1 ttt2]);
            TT(i,:) = fliplr(-ttt);
            K1(i,:)=G1(ttt,kv);
            K2(i,:)=G2(ttt,kv);
        else
            ttt1 = logspace(log10(alpmin*(((kappa*epsi)/(kappac))^2/c)),log10(alp*(((kappa*epsi)/(kappac))^2/c)),KL/2);
            ttt2 = linspace((alpmin*(((kappa*epsi)/(kappac))^2/c)),(alp*(((kappa*epsi)/(kappac))^2/c)),KL/2+2);
            ttt2 = ttt2(2:end-1);
            ttt = sort([ttt1 ttt2]);
            TT(i,:) = fliplr(-ttt);
            K1(i,:)= 0.; % G1(ttt,kv+eps);
            K2(i,:)= -c.*kappac./kappa./epsi./sqrt(pi.*c.*ttt) + ...
                    c.*(kappac./kappa./epsi).^2*exp((kappac./kappa./epsi)^2*c.*ttt) ...
                    .*erfc(kappac./kappa./epsi.*sqrt(c.*ttt));
            % (- c*k.^2.* (1 + FF(k)).*( exp(-c*k.^2*t)./sqrt(pi*c*k.^2*t) - FF(k).*exp((FF(k).^2 - 1).*c*k.^2*t).*erfc(FF(k).*sqrt(c*k.^2*t)) )       );
% G2(ttt,kv+eps);
        end



    end
    K1(isnan(K1) | isinf(K1)) = 0;
    K2(isnan(K2) | isinf(K2)) = 0;

    K1 = fliplr(K1);
    K2 = fliplr(K2);

    %%
    Fs = 1./(x(2)-x(1));
    fre = Fs*(-N/2:1:N/2-1)/N;
    kv = 2*pi*fre'; %!

    filter = (sin(pi*(abs([-N/2:1:N/2-1]))/(N*0.5))./(pi*(abs([-N/2:1:N/2-1]))/(N*0.5))).^(1)';
    filter(isnan(filter)) = 1;

    % F/(F+1)
    FFoFF = FF(kv+eps)./(FF(kv+eps)+1);

    FFoFF(isinf(FFoFF) | isnan(FFoFF)) = 1;

    if Elastic_Flag == 1
        FFoFF(:) = 1.;
    end
    %% various things are computed outside the loop to save time:
    dhat =     exp(1i*kv.*(x-0.5*((N+1)/Fs))); %!
    dhat2 =    filter.*exp(1i*(kv.*(x-0.5*((N+1)/Fs))));
    akdhat1=   -N*G/(2*(1-nuu))*abs(kv);
    akdhat2 =  -N*1i*G*B/3.*(1+nuu)/(1-nuu).*kv.*filter;
    akdhat3 =  -N*FFoFF;
    akdhat4 =   N*3./(2*B*(1+nuu)).*FFoFF;
    akdhat5 =   N*G*B/3.*(1+nuu)/(1-nuu).*abs(kv);

    % For parallel purposes
    TT_K = [TT(:,end).*K1(:,end), TT(:,end).*K2(:,end)];
    K_var(:,:,1) = K1; K_var(:,:,2) = K2;
    TT_diff = TT(:,2:end) - TT(:,1:end-1);
    con_var = zeros(size(dx,1),6);
    DD_var = zeros(size(dx,1), KL, 3);
    ind_0 = [1,2,2,2,1,1];
    ind_1 = [1,1,2,3,3,2];
    ind_2 = [1,2,2,2,1,1];
    ind_3 = [1,1,2,3,3,2];
    ind_4 = [1,2,2,2,1,1];
    ind_5 = [1,1,2,3,3,2];
    %% CHECK SEPARATION OF SCALES
    if epsi*2*pi/(x(2)-x(1)) > 0.05
        error('epsi is too large compered to the mininum resolved wavelength: make epsi larger or the gridsize larger')
    end
    %% Main Loop Starts
    tryagain = 0;

    Lb = G/(1-nuu)*L/(b(1)*si0);

    Linf = Lb/pi*(b(1)/(b(1)-a(1)))^2;

    counter = 0;
    maxtolvio = 0;
    cumutolvio = 0;

    % Main loop cannot be parallelized
    for it = 2:NT
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
            dphi = dphi0 - gamma*log(Vr.*theta/L);
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
            
            % IF elastic, no cross-fault fluid motion is allowed
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
                            if gamma ~= 0
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
                    filename = strcat('../outputMats/', 'Elastic_Flag', num2str(Elastic_Flag), '_FluxTime_', num2str(flux), '_NewFH_', num2str(FHFlag), '_nu_nuu_', num2str(nu), '_',  num2str(nuu), '_gamma_', num2str(gamma),...
                                      '_pflag_', num2str(poreflag),'_c_', num2str(cc), '_factors_', ...
                                      num2str(factors(1)), '_', num2str(factors(2)), '_',num2str(factors(3)), '_', num2str(Terminating_time), '_', num2str(initial_V_ratio), '_312.mat');

                    % Record excuting time of the program
                    %t1 = cputime - t0;
                    t1 = toc(tstart);
                    disp('Finished!');
                    save(filename);

                    % Write changable parameters into a '.txt' file
                    txtname = strcat('../outputMats/', 'Elastic_Flag', num2str(Elastic_Flag), '_FluxTime_', num2str(flux), '_NewFH_', num2str(FHFlag), '_nu_nuu_', num2str(nu), '_',  num2str(nuu), '_gamma_', num2str(gamma),...
                                      '_pflag_', num2str(poreflag),'_c_', num2str(cc), '_factors_', ...
                                      num2str(factors(1)), '_', num2str(factors(2)), '_',num2str(factors(3)),'_', num2str(Terminating_time), '_', num2str(initial_V_ratio), '.txt');


                    fileID = fopen(txtname, 'w');

                    fprintf(fileID, '%25s', 'Filename: '); fprintf(fileID, filename);
                    fprintf(fileID, '\n%25s', 'NT: '); fprintf(fileID, num2str(NT));
                    fprintf(fileID, '\n%25s', 'lhs, rhs: '); 
                    fprintf(fileID, strcat(num2str(lhs), ', ', num2str(rhs))); 
                    fprintf(fileID, '\n%25s', 'Gamma: '); fprintf(fileID, num2str(gamma));
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
                    disp(strcat('Gamma = ', num2str(gamma)));
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
    filename = strcat('../outputMats/', 'Elastic_Flag', num2str(Elastic_Flag), '_FluxTime_', num2str(flux), '_NewFH_', num2str(FHFlag), '_nu_nuu_', num2str(nu), '_',  num2str(nuu), '_gamma_', num2str(gamma),...
                      '_pflag_', num2str(poreflag),'_c_', num2str(cc), '_factors_', ...
                      num2str(factors(1)), '_', num2str(factors(2)), '_',num2str(factors(3)),'_', num2str(Terminating_time), '_', num2str(initial_V_ratio), '_312.mat');


    % Record excuting time of the program
    % t1 = cputime - t0;
    t1 = toc(tstart);
    disp('Finished');
    save(filename);
    
    % Write changable parameters into a '.txt' file
    txtname = strcat('../outputMats/', 'Elastic_Flag', num2str(Elastic_Flag), '_FluxTime_', num2str(flux), '_NewFH_', num2str(FHFlag), '_nu_nuu_', num2str(nu), '_',  num2str(nuu), '_gamma_', num2str(gamma),...
                      '_pflag_', num2str(poreflag),'_c_', num2str(cc), '_factors_', ...
                      num2str(factors(1)), '_', num2str(factors(2)), '_',num2str(factors(3)),'_', num2str(Terminating_time), '_', num2str(initial_V_ratio), '.txt');

    
    fileID = fopen(txtname, 'w');
    
    fprintf(fileID, '%25s', 'Filename: '); fprintf(fileID, filename);
    fprintf(fileID, '\n%25s', 'NT: '); fprintf(fileID, num2str(NT));
    fprintf(fileID, '\n%25s', 'lhs, rhs: '); 
    fprintf(fileID, strcat(num2str(lhs), ', ', num2str(rhs))); 
    fprintf(fileID, '\n%25s', 'Gamma: '); fprintf(fileID, num2str(gamma));
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

 % 1-0-1-0 injection, period T = 2020 s. 
function InMass = INjectmass(t)
    flux = 1.e-4;
    if ((t / 2020.) - floor(t / 2020.)) <= 0.5
        InMass = flux * 0.5 * floor(t / 2020.) * 2020. + ... 
                 flux * ((t / 2020.) - floor(t / 2020.)) * 2020.; 
    else
        InMass = flux * 0.5 * ceil(t / 2020.) * 2020.;
    end
end
