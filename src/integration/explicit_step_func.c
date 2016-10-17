void explicit_step_func(double aplanet, double *y,double *mdL, double *mdR) {
    int i;
    
    double TL = 0;
    double TR = 0;
    double facm,facp,am,ap,bm,bp,nup,num,drm,drp;
    double fac,rm,rp,mfacL,mfacR,lfacL,ufacR;
    double yL,yM,yR;

    double rgL = .5*(rmin[0] + exp(log(rmin[0]) - dlr));
    double rgR = .5*(rmin[NR] + exp(log(rmin[NR]) + dlr));
    double bc_fac_L = 2*M_PI*rgL;
    double bc_fac_R = 2*M_PI*rgR;

    for(i=0;i<NR;i++) {

        fac = dr[i]*dTr_ex(rc[i],aplanet)*y[i];
        if (rc[i] < aplanet) TL += fac;
        else    TR += fac;
    }

    for(i=0;i<NR;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        facm = dep_func(rm,aplanet,planet.xd,planet.wd);
        facp = dep_func(rp,aplanet,planet.xd,planet.wd);
        facm *= (rm < aplanet) ? TL : TR; 
        facp *= (rp > aplanet) ? TR : TL; 

        num = nu(rm);
        nup = nu(rp);

        am = 3*num;
        bm = 3*(num/rm)*(params.gamma - .5);

        ap = 3*nup;
        bp = 3*(nup/rp)*(params.gamma-.5);

        drm = rc[i]-rc[i-1];
        drp = rc[i+1]-rc[i];

        mfacR = ( bp * ( rc[i+1] - rp) - ap) / drp;
        mfacL = ( am + bm * ( rm - rc[i-1] ) ) / drm;
        ufacR  = ( ap + bp * ( rp - rc[i] ) ) / drp;
        lfacL  = ( bm *(rc[i] - rm ) - am ) / drm; 
        yL = (i==0) ? params.bc_val[0]*bc_fac_L: y[i-1] ;

        yM = y[i];
        yR = (i==NR-1) ? params.bc_val[2]*bc_fac_R: y[i+1];
            
        mdL[i] = (lfacL*yL + mfacL*yM) - 2*sqrt(rp)*facp;
        mdR[i] = (mfacR*yM + ufacR*yR) - 2*sqrt(rm)*facm;
        //dy[i] = dt*(lfac*yL + mfac*yM + ufac*yR) ;
        //dy[i] += dt*( -2*sqrt(rm)*facm + 2*sqrt(rp)*facp);


    }

    return;
}
