# MUN SARS-CoV-2 modelling team
# Branching random walk model of infection

# Parameters of the distributions considered
# Proportion of subclinical infections: probabilita.subclinica
# Proportion of severe infection given clinical condition: cond.prob.accuta
# Generation time distribution: (Weibull) tempo.di.generazione
# From infection to recover subclinical: (Gamma) tempo.di.recuperare.subcli
# From infection to onset: (Gamma) tempo.di.insorgenza
# From onset to recover non-accute: (Gamma) tempo.di.recuperare.nonaccuta
# From onset to hospitalisation: (Exponential)  tasso.ammisione
# From hospital admission to dimission: (Exponential) tasso.dimissione
# From Hospital dimission to recover: (Exponential) tasso.di.recupero
# Infected importation rate: tasso.di.immigrazione.infetto
para.infezione<-list(probabilita.subclinica=0.15,
                     tempo.di.generazione=c(forma=2.83,scala=5.67),
                     tasso.recuperare.subcli=0.06,
                     tasso.recuperare.nonaccuta=0.08,
                     tempo.di.insorgenza=c(forma=5.8,scala=0.95),
                     cond.prob.accuta=0.078,
                     tasso.ammisione=0.35,
                     tasso.dimissione=0.1,
                     tasso.di.recupero=0.25,
                     tasso.di.immigrazione.infetto=1/30
                     )

# inefectivity
# Relative infectiousness of preclinical, clinical and subclinical 
# conditions, in this case with respect to preclinical (fully ineffective)
infettivita<-list(preclinica=1.0, clinica=0.3, subclinica=0.75)
# intervention as a function of time (day), size of the jump and
#  changepoint time.
intervento<-function(tempo, dimensione.del.salto, punto.di.cambiamento){
    for(i in seq(length(punto.di.cambiamento),1))
        if(tempo >= punto.di.cambiamento[i]) return(dimensione.del.salto[i])
}

# new infections for the day
nuove.infezioni<-function(oggi,infettatore,infetti,distribuzioni){
    caso<-data.frame(stringsAsFactors = FALSE)
    n0<-length(infettatore)
    stopifnot(n0 == length(infetti))
    for(i in seq(1,n0)){
        if(infetti[i] <1) next
        for(j in seq(1,infetti[i])){
            if(distribuzioni$probabilita.subclinica > runif(1,0,1)){
                dor<-oggi+rexp(1,
                       rate=distribuzioni$tasso.recuperare.subcli)
                caso<-rbind(caso,c(oggi,NA,NA,NA,dor,infettatore[i]))
            }else{
                doo<-oggi+rgamma(1,
                           shape=distribuzioni$tempo.di.insorgenza['forma'],
                           scale=distribuzioni$tempo.di.insorgenza['scala'])
                if(distribuzioni$cond.prob.accuta > runif(1,0,1)){
                        doa<-doo+rexp(1,rate=distribuzioni$tasso.ammisione)
                        dod<-doa+rexp(1,rate=distribuzioni$tasso.dimissione)
                        dor<-dod+rexp(1,rate=distribuzioni$tasso.di.recupero)
                        caso<-rbind(caso,c(oggi,doo,doa,dod,dor,infettatore[i]))
                }
                else{
                  rec<-rexp(1,
                     rate=distribuzioni$tasso.recuperare.nonaccuta)
                  caso<-rbind(caso,c(oggi,doo,NA,NA,doo+rec,infettatore[i]))
                }
            }
        }
    }
    if(dim(caso)[1]>0) names(caso)<-c('doi','doo','doa','dod','dor','infector')
    return(caso)
}



# for a given day, split the data into preclinical clinical and subclinical
# and call the nuove.infezioni() to generate the infections of the day
infezione_oggi<-function(oggi,casi,lambda,para.infezione,infettivita,omega.gt){
    rownames(casi)<-seq(1,nrow(casi))
    names(casi)<-c('doi','doo','doa','dod','dor','infector')   
    pgio<-nuove.infezioni(oggi,0,
                          infetti=rpois(1,
                            para.infezione$tasso.di.immigrazione.infetto),
                          para.infezione)
    precli<-casi[!is.na(casi$doo) & casi$doi>(oggi-21) & casi$doo > oggi,]
    if(dim(precli)[1]>0){
        lambda0<-omega.gt[ceiling(oggi-precli$doi)]*
            lambda*infettivita$preclinica
         npre=nuove.infezioni(oggi,infettatore=as.numeric(rownames(precli)),
                             infetti=rpois(length(lambda0),lambda0),
                             para.infezione)
     } else npre<-data.frame()
    clinic<-casi[!is.na(casi$doa) & casi$doi>(oggi-21) &
                 casi$doa <= oggi & casi$dor >= oggi,]
    if(dim(clinic)[1]>0){
        lambda0<-omega.gt[ceiling(oggi-clinic$doi)]*
            lambda*infettivita$clinica
         ncli=nuove.infezioni(oggi,infettatore=as.numeric(rownames(clinic)),
                             infetti=rpois(length(lambda0),lambda0),
                             para.infezione)
    }else ncli<-data.frame()
    subcli<-casi[is.na(casi$doo) & casi$doi>(oggi-21) & casi$dor >= oggi,]
    if(dim(subcli)[1]>0){
        lambda0<-omega.gt[ceiling(oggi-subcli$doi)]*
            lambda*infettivita$subclinica
        nsub=nuove.infezioni(oggi,infettatore=as.numeric(rownames(subcli)),
                             infetti=rpois(length(lambda0),lambda0),
                             para.infezione)
    } else nsub<-data.frame()
    pgio<-rbind(pgio,npre,ncli,nsub)
    return(pgio)
}

# main function for the branching random walk
camminata.casuale<-function(base.R0, punto.di.cambia, i0,
                         para.infezione,infettivita,
                         poppolazione=5.5e5, giorni.totali=210){
    controllo<-base.R0/base.R0[1]
    Ro<-base.R0[1]/(1-(1-infettivita$subclinica)
        * para.infezione$probabilita.subclinica)
    casi<-nuove.infezioni(0,0,i0,para.infezione)
    names(casi)<-c('doi','doo','doa','dod','dor','infector')       
    omega.gt<-pweibull(seq(1, 21, 1),
                shape = para.infezione$tempo.di.generazione['forma'],
                scale = para.infezione$tempo.di.generazione['scala'])
    omega.gt<-omega.gt-c(0,omega.gt[-21])
    omega.gt<-omega.gt/sum(omega.gt)
    for(oggi in seq(1,giorni.totali)){
        n0.casi<-nrow(casi)
        nuovi.casi<-infezione_oggi(oggi,casi,(1-n0.casi/poppolazione)*
                        Ro*intervento(oggi,controllo,punto.di.cambia),
                        para.infezione,infettivita,omega.gt)
         casi<-rbind(casi,nuovi.casi)
    }
    return(casi)
}


per_giorno<-function(casi, giorni.totali=210){
    names(casi)<-c('doi','doo','doa','dod','dor','infector')    
    modello<-matrix(0,nrow=giorni.totali+1,ncol=8)
    # infected (deltaSI)
    modello[,1]<-hist(casi$doi,breaks = seq(-1,giorni.totali),
                      plot = FALSE)$count
    # preclinical (deltaSP)
    modello[,2]<-hist(casi$doi[!is.na(casi$doo)],breaks = seq(-1,giorni.totali),
                      plot = FALSE)$count
    # subclinical (deltaSS)
    modello[,3]<-hist(casi$doi[is.na(casi$doo)],breaks = seq(-1,giorni.totali),
                      plot = FALSE)$count
    # clinical (deltaPC)
    modello[,4]<-hist(casi$doo[!is.na(casi$doo) & casi$doo <= giorni.totali],
                      breaks = seq(-1,giorni.totali), plot = FALSE)$count
    # admitted to hospital (deltaCH)
    modello[,5]<-hist(casi$doa[!is.na(casi$doa) & casi$doa <= giorni.totali],
                      breaks = seq(-1,giorni.totali), plot = FALSE)$count
    # dismissed from hospital (deltaHC)
    modello[,6]<-hist(casi$dod[!is.na(casi$doa) & casi$dod <= giorni.totali],
                      breaks = seq(-1,giorni.totali), plot = FALSE)$count
    # clinical recovered (deltaCR) collects both 'accute' and 'nonaccute'
    modello[,7]<-hist(casi$dor[!is.na(casi$doo) & casi$dor <= giorni.totali],
                      breaks = seq(-1,giorni.totali), plot = FALSE)$count    
    # subclinical recovered (deltaSR)
    modello[,8]<-hist(casi$dor[is.na(casi$doo) & casi$dor <= giorni.totali],
                      breaks = seq(-1,giorni.totali), plot = FALSE)$count
    modello<-as.data.frame(modello)
    names(modello)<-c('deltaSI','deltaSP','deltaSS','deltaPC',
                      'deltaCH','deltaHC','deltaCR','deltaSR')
    return(modello)
}


branching.3ihr<-function(cases.per.day, total.days){
    sir<-data.frame(C=cumsum(cases.per.day$deltaSI),
                    I=cumsum(cases.per.day$deltaSI-cases.per.day$deltaCR
                             - cases.per.day$deltaSR),
                    Is=cumsum(cases.per.day$deltaSS-cases.per.day$deltaSR),
                    Ip=cumsum(cases.per.day$deltaSP-cases.per.day$deltaPC),
                    Ic=cumsum(cases.per.day$deltaPC-cases.per.day$deltaCR),
                    H=cumsum(cases.per.day$deltaCH-cases.per.day$deltaHC),
                    R=cumsum(cases.per.day$deltaCR+cases.per.day$deltaSR))
    return(sir)
}


facciamolo<-function(R0,changepoint,i0,
                     para.infezione,infettivita,
                     population.size=5.5e5, Days=210,
                     Runs=100,ymax=150){

    time<-seq(0,Days)
    plot(c(0,Days),c(0,ymax),col='blue',xlab='day',
         ylab='count',type='n',xlim=c(0,Days),ylim=c(0,ymax))
    oi1<-rep(0,Days)
    oc1<-rep(0,Days)
    os1<-rep(0,Days)
    op1<-rep(0,Days)
    oh1<-rep(0,Days)
    or1<-rep(0,Days)
    for(i in 1:Runs){
        cases<-camminata.casuale(base.R0=R0,
                         punto.di.cambia=changepoint, i0=i0,
                         para.infezione=para.infezione,
                         infettivita=infettivita,
                         poppolazione=population.size, giorni.totali=Days)
        cases.per.day<-per_giorno(cases, giorni.totali=Days)
        risultato<-branching.3ihr(cases.per.day, total.days=Days)
        
        lines(time,risultato$C,col='pink')
        lines(time,risultato$Ip,col='green')
        lines(time,risultato$Ic,col='blue1')
        lines(time,risultato$Is,col='grey')
        lines(time,risultato$H,col='red')
        #if (!exists(ocl)) {
        #    
        #}
        oc1<-oc1+risultato$C
        op1<-op1+risultato$Ip
        oi1<-oi1+risultato$Ic
        os1<-os1+risultato$Is
        oh1<-oh1+risultato$H
        #print(length(risultato$Rs))
        if (length(risultato$Rs) != 0) {
            or1<-or1+risultato$Rs
        }
        #or1<-or1+risultato$Rs FIX THIS
    }
    or1<-(oc1-or1-os1)/Runs
    oi1<-oi1/Runs
    oc1<-oc1/Runs
    os1<-os1/Runs
    oh1<-oh1/Runs
    op1<-op1/Runs

    lines(time,oc1,col='magenta',lwd=2)
    lines(time,op1,col='green4',lwd=2)
    lines(time,oi1,col='blue4',lwd=2)
    lines(time,os1,col='grey4',lwd=2)
    lines(time,oh1,col='red',lwd=2)

#lines(1:length(nl0),nl0,col='brown')
#points(1:length(nl0),nl0)
#lines(1:length(nl0),active0,)
    lines(time,or1, lwd=2) #fix????
}


#=================================================================
#
cases<-camminata.casuale(base.R0=c(3.6,0.5,1.6,0.5),
                         punto.di.cambia=c(1,15,60,120), i0=4,
                         para.infezione,infettivita,
                         poppolazione=5.5e5, giorni.totali=210)
cases.per.day<-per_giorno(cases, giorni.totali=210)
branching.3ihr(cases.per.day, total.days=210)
#C:\Users\colem\Desktop\Infectious Disease Research\JC_RModel
png(filename='C:/Users/colem/Desktop/Infectious Disease Research/JC_RModel/terranova0.png',width=800,height=640)
facciamolo(R0=c(4.6,0.5,1.6,0.5),changepoint=c(1,15,60,120),i0=4,
           para.infezione,infettivita,
           population.size=5.5e5, Days=210,
           Runs=100,ymax=150)
legend(120,120,c('cumulative','preclinical','clinical','subclinical','hospitalised','reported'),lty=1,lwd=2,col=c('magenta','green4','blue4','grey4','red','black'))  
dev.off()



