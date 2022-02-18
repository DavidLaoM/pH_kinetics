% % appendix_ODEtranslate2latex

%% PFK
clc

syms Vm
v_PFK =  Vm;
latex(v_PFK)

syms Vm FBP GAP DHAP Keq
v_ALD =  Vm .* (FBP-(GAP.*DHAP)./ Keq);
latex(v_ALD)

syms Vm DHAP NADH G3P NAD Keq
v_GPD = Vm .* (DHAP .* NADH - G3P .* NAD ./ Keq);
latex(v_GPD)

syms Vm DHAP GAP Keq
v_TPI = Vm .* (DHAP - GAP ./ Keq);
latex(v_TPI)

syms Vm gR F6P Kf6p ATP Katp L Ciatp Kiatp Kiatp Camp AMP Kamp Cf26bp F26BP Kf26bp Cf16bp F16BP Kf16bp Catp
PFK_nom=(Vm.*gR.*(F6P./Kf6p).*(ATP./Katp).*(1+(F6P./Kf6p)+(ATP./Katp)+gR.*((F6P./Kf6p).*(ATP./Katp))));
PFK_denom=(1+F6P./Kf6p+ATP./Katp+(gR.*(F6P./Kf6p).*(ATP./Katp))).^2+...
    L.*...
    ((1+Ciatp.*(ATP./Kiatp))./(1+ATP./Kiatp)).^2.*...
    ((1+Camp.*(AMP./Kamp))./(1+AMP./Kamp)).^2.*...
    ((1+((Cf26bp*F26BP)./(Kf26bp))+((Cf16bp.*F16BP)./(Kf16bp)))./(1+(F26BP./Kf26bp)+(F16BP./Kf16bp))).^2.*...
    (1+Catp.*(ATP./Katp)).^2;
v_PFK_long = PFK_nom ./ PFK_denom;
latex(v_PFK_long)


%% HXK
clc

syms Vm ATP GLCi ADP G6P Keq Katp Kglc Kadp Kg6p
v_GLK = (Vm .* (ATP.*GLCi-((ADP.*G6P)./Keq))) ./ ...
    ((Katp.*Kglc) .* (1+ATP./Katp+ADP./Kadp) .* ...
    (1+GLCi./Kglc+G6P./Kg6p));
latex(v_GLK)

syms Vm NADP G6P NADPH PG6 Keq
v_G6PDH = Vm .* (NADP .* G6P - NADPH .* PG6 ./ Keq);
latex(v_G6PDH)

%% GAPDHliterature
clc

syms BPG NADH NAD GAP Kbpg Knadh Kgap Knad Vmf Vmr
v_GAPDH_lit = (-(Vmr .* BPG .* NADH ./ (Kbpg .* Knadh)) + Vmf .* NAD .* GAP ./ ( Kgap .* Knad)) ./ ((1 + NAD ./ Knad + NADH ./ Knadh) .* (1 + BPG ./ Kbpg + GAP ./ Kgap));
latex(v_GAPDH_lit)


%% GAPDHfwd
clc

syms Vm BPG ADP P3G ATP Keq
v_PGK(Vm,BPG,ADP,P3G,ATP,Keq) = Vm .* (BPG .* ADP - P3G .* ATP ./ Keq);
latex(v_PGK)

syms Vmf BPG NADH GAP NAD Keq Kbpg Knadh Kgap Knad
v_GAPDHfwd = (Vmf .* (GAP .* NAD - BPG .* NADH ./ Keq)./(Kgap .* Knad))./...
    ((1 + NAD ./ Knad + NADH ./ Knadh) .* (1 + BPG ./...
    Kbpg + GAP ./ Kgap));
latex(v_GAPDHfwd)


%% GAPDHrev
clc

syms Vm BPG ADP P3G ATP Keq
v_PGK(Vm,BPG,ADP,P3G,ATP,Keq) = Vm .* (BPG .* ADP - P3G .* ATP ./ Keq);
latex(v_PGK)

syms Vmr BPG NADH GAP NAD Keq Kbpg Knadh Kgap Knad
v_GAPDHrev = ((Vmr .* (BPG .* NADH ./ Keq - GAP .* NAD))./(Kbpg .* Knadh))./...
    ((1 + NAD ./ Knad + NADH ./ Knadh) .* (1 + BPG ./...
    Kbpg + GAP ./ Kgap));
latex(v_GAPDHrev)


%% ENO
clc

syms P2G PEP K2pg Keq Kpep vm vENO
v_ENO(P2G,PEP,K2pg,Keq,Kpep,vm)=((vm./K2pg).*(P2G-PEP./Keq))./...
    (1+P2G./K2pg+PEP./Kpep);
latex(v_ENO)


%% PYK
clc

syms Vm Kadp Kpep  ADP PEP  ADP Kadp PEP Kpep PEP Kpep   hill L   ATP Katp    FBP Kf16bp  hill  PEP Kpep   hill
v_PYK(Vm,Kadp,Kpep,ADP,PEP,hill,L,ATP,Katp,FBP,Kf16bp)=(((Vm./(Kadp.*Kpep).*ADP.*PEP)./((1+ADP./Kadp).*(1+PEP./Kpep))).*((PEP./Kpep+1).^hill./(L.*((ATP./Katp+1)./(FBP./Kf16bp+1)).^    hill+(PEP./Kpep+1).^hill)));
latex(v_ENO)

syms Vm PYR NADH LAC NAD Keq
v_LDH(Vm,PYR,NADH,LAC,NAD,Keq) = Vm .*(PYR .* NADH - LAC .* NAD ./ Keq);
latex(v_LDH)


%% PDC
clc

syms vm PYR Kpyr hill
v_PDC(vm,PYR,Kpyr,hill) = (vm.*(PYR./Kpyr).^hill)./...
        (1+(PYR./Kpyr).^hill);
latex(v_PDC)

syms Vm AcAld NADH ETOH NAD Keq
v_ADH(Vm,AcAld,NADH,ETOH,NAD,Keq) = Vm .* (AcAld .* NADH - ETOH .* NAD ./ Keq);
latex(v_ADH)

