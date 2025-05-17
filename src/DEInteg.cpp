#include "../include/DEInteg.h"
#include "../include/SAT_Const.h"
#include "../include/sign_.h"
#include <cmath>
#include <iostream>
#include <iomanip>
    /**
     * @file DEInteg.cpp
     * @brief El archivo contiene las implementaciones de DEInteg.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
Matrix& DEInteg(Matrix& f(double t,Matrix z),double t, double tout,double relerr,double abserr,int n_eqn,Matrix &y){
  if(y.n_row<y.n_column){
    y=transpose(y);
  }
  cout<<setprecision(15);

  Matrix yout=zeros(n_eqn,1),ypout=zeros(n_eqn,1),two(14),gstr(14);
  bool start=false,phase1=false,nornd=false,crash=false,success=false,PermitTOUT=false,OldPermit=false,stiff=false;
  long double delsgn=0.0,x=0.0,hi=0.0,ki=0.0,kold=0.0,temp1=0.0,term=0.0,psijm1=0.0,eta=0.0,sum=0.0,absh=0.0,hold=0.0,hnew=0.0,
  k=0.0,round=0.0,gamma=0.0,i=0.0,p5eps=0.0,ifail=0.0,kp1=0.0,kp2=0.0,km1=0.0,km2=0.0,ns=0.0,nsp1=0.0,realns=0.0,im1=0.0,temp2=0.0,
  temp3=0.0,reali=0.0,temp4=0.0,nsm2=0.0,limit1=0.0,temp5=0.0,temp6=0.0,limit2=0.0,nsp2=0.0,ip1=0.0,tau=0.0,xold=0.0,erkm2=0.0,
  erkm1=0.0,erk,err=0.0,knew=0.0,rhi=0.0,h=0.0,erkp1=0.0,rhodouble=0.0,told=0.0,epsilon=0.0,del=0.0,absdel=0.0,tend=0.0,nostep=0.0,kle4=0.0,
  releps=0.0,abseps=0.0,twou=0.0,fouru=0.0;
  double r=0.0;
  int l=0;
  twou  = 2*eps;
  fouru = 4*eps;


  struct DE_STATE_t {
    int DE_INIT = 1;      
    int DE_DONE = 2;      
    int DE_BADACC = 3;    
    int DE_NUMSTEPS = 4;  
    int DE_STIFF = 5;     
    int DE_INVPARAM = 6;  
  };

  DE_STATE_t DE_STATE;

  int State_ = DE_STATE.DE_INIT;
PermitTOUT = true,OldPermit;         
told = 0;


double arrtwo[]  = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};
double arrgstr[] = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679,0.00592, 0.00524, 0.00468};
for(int wxy=1;wxy<=14;wxy++){
  two(wxy)=arrtwo[wxy-1];
  gstr(wxy)=arrgstr[wxy-1];
}

Matrix yy    = zeros(n_eqn,1);    
Matrix wt    = zeros(n_eqn,1);
Matrix p     = zeros(n_eqn,1);
Matrix yp    = zeros(n_eqn,1);
Matrix phi   = zeros(n_eqn,17);
Matrix g     = zeros(14,1);
Matrix sig   = zeros(14,1);
Matrix rho   = zeros(14,1);
Matrix w     = zeros(13,1);
Matrix alpha = zeros(13,1);
Matrix beta  = zeros(13,1);
Matrix v     = zeros(13,1);
Matrix psi_  = zeros(13,1);

if (t==tout) {   
  return y;}

  epsilon = fmax(relerr,abserr);
  if ( ( relerr <  0.0  ) || ( abserr <  0.0 ) ||  ( epsilon    <= 0.0  ) ||
    ( State_  >  DE_STATE.DE_INVPARAM ) ||  ( (State_ != DE_STATE.DE_INIT) &&  (t != told)           ) )
  {
     State_ = DE_STATE.DE_INVPARAM;  
	   return y;                                  
  }

  del    = tout - t;
  absdel = fabs(del);

  tend   = t + 100.0*del;
  if (!PermitTOUT){
    tend = tout;
  }
  nostep = 0;
  kle4   = 0;
  stiff  = false;
  releps = relerr/epsilon;
  abseps = abserr/epsilon;

  if  ( (State_==DE_STATE.DE_INIT) || (!OldPermit) || (delsgn*del<=0.0) ){
    
    
    start  = true;
    x      = t;
    yy     = y;
    delsgn = sign_(1.0, del);
    h = sign_( fmax(fouru*fabs(x), fabs(tout-x)), tout-x );}
while(true){
  if (fabs(x-t) >= absdel){
    yout  = zeros(n_eqn,1);
    ypout = zeros(n_eqn,1);
    g(2)   = 1.0;
    rho(2) = 1.0;
    hi = tout - x;
    ki = kold + 1;

    for (int i=1;i<=k;i++){
      temp1 = i;
      w(i+1) = 1.0/temp1;}
      
      term = 0.0;
      for (int j=2;j<=ki;j++){
        psijm1 = psi_(j);
        gamma = (hi + term)/psijm1;
        eta = hi/psijm1;
        for (int i=1;i<=ki+1-j;i++){
          w(i+1) = gamma*w(i+1) - eta*w(i+2);}
          g(j+1) = w(2);
          rho(j+1) = gamma*rho(j);
          term = psijm1;}
      
        if(yout.n_row>yout.n_column){
        yout=transpose(yout);}
        if(ypout.n_row>ypout.n_column){
        ypout=transpose(ypout);}
        if(y.n_row>y.n_column){
        y=transpose(y);
         }
          for (int j=1;j<=ki;j++){
            i = ki+1-j;
            yout  = yout  + extract_column(phi,i+1)*g(i+1);
            ypout = ypout + extract_column(phi,i+1)*rho(i+1);
          }
            yout = y +yout*hi;
            y    = yout;
      State_    = DE_STATE.DE_DONE; 
      t         = tout;             
      told      = t;                
      OldPermit = PermitTOUT;
      return y;                       
    }          

    if ( !PermitTOUT && ( fabs(tout-x) < fouru*fabs(x) ) ){
      h = tout - x;
      yp = f(x,yy);          
      y = yy + yp*h;                
      State_    = DE_STATE.DE_DONE; 
      t         = tout;             
      told      = t;                
      OldPermit = PermitTOUT;
      return y;                       
    }

    h  = sign_(min(fabs(h), fabs(tend-x)), h);
    for (l=1;l<=n_eqn;l++){
      wt(l) = releps*fabs(yy(l)) + abseps;
    }

    if (fabs(h) < fouru*fabs(x)){
      h = sign_(fouru*fabs(x),h);
      crash = true;
      return y;
  }

  p5eps  = 0.5*epsilon;
  crash  = false;
  g(2)   = 1.0;
  g(3)   = 0.5;
  sig(2) = 1.0;

  ifail = 0;

  round = 0.0;
  for (l=1;l<=n_eqn;l++){
    round = round + (y(l)*y(l))/(wt(l)*wt(l));
  }
  round = twou*sqrt(round);
  if (p5eps<round){
    epsilon = 2.0*round*(1.0+fouru);
    crash = true;
    return y;
  }
  if (start){
  
    yp = transpose(f(x,y));
    sum = 0.0;
    for (l=1;l<=n_eqn;l++){
      phi(l,2) = yp(l);
      phi(l,3) = 0.0;
      sum = sum + (yp(l)*yp(l))/(wt(l)*wt(l));
    }
    sum  = sqrt(sum);
    absh = fabs(h);
    if (epsilon<16.0*sum*h*h){
      absh=0.25*sqrt(epsilon/sum);
    }
    h    = sign_(fmax(absh, fouru*fabs(x)), h);
    hold = 0.0;
    hnew = 0.0;
    k    = 1;
    kold = 0;
    start  = false;
    phase1 = true;
    nornd  = true;
    if (p5eps<=100.0*round){
      nornd = false;
      for (l=1;l<=n_eqn;l++){
        phi(l,16)=0.0;
      }
    }
  }

  while(true){

    kp1 = k+1;
    kp2 = k+2;
    km1 = k-1;
    km2 = k-2;

    if (h !=hold){
      ns=0;
    }
    if (ns<=kold){
      ns=ns+1;
    }
    nsp1 = ns+1;

    if (k>=ns){
      beta(ns+1) = 1.0;
      realns = ns;
      alpha(ns+1) = 1.0/realns;
      temp1 = h*realns;
      sig(nsp1+1) = 1.0;
      if (k>=nsp1){
        for (int i=nsp1;i<=k;i++){
          im1   = i-1;
          temp2 = psi_(im1+1);
          psi_(im1+1) = temp1;
          beta(i+1)  = beta(im1+1)*psi_(im1+1)/temp2;
          temp1    = temp2 + h;
          alpha(i+1) = h/temp1;
          reali = i;
          sig(i+2) = reali*alpha(i+1)*sig(i+1);
        }
      }
      psi_(k+1) = temp1;
      
      
      if (ns>1){
          
        if (k>kold){
          temp4 = k*kp1;
          v(k+1) = 1.0/temp4;
          nsm2 = ns-2;
          for (int j=1;j<=nsm2;j++){
            i = k-j;
            v(i+1) = v(i+1) - alpha(j+2)*v(i+2);
          }
        }

          
        limit1 = kp1 - ns;
        temp5  = alpha(ns+1);
        for (int iq=1;iq<=limit1;iq++){
          v(iq+1) = v(iq+1) - temp5*v(iq+2);
          w(iq+1) = v(iq+1);
        }
        g(nsp1+1) = w(2);}
        else{
          for (int iq=1;iq<=k;iq++){
            temp3 = iq*(iq+1);
            v(iq+1) = 1.0/temp3;
            w(iq+1) = v(iq+1);
          }
        }

      
        nsp2 = ns + 2;
        if (kp1>=nsp2){
          for (int i=nsp2;i<=kp1;i++){
            limit2 = kp2 - i;
            temp6  = alpha(i);
            for (int iq=1;iq<=limit2;iq++){
              w(iq+1) = w(iq+1) - temp6*w(iq+2);
            }
            g(i+1) = w(2);
          }
        }
  } 
  
  
  if (k>=nsp1){
    for (int i=nsp1;i<=k;i++){
      temp1 = beta(i+1);
      for (l=1;l<=n_eqn;l++){
        phi(l,i+1) = temp1 * phi(l,i+1);
      }
    }
  }
  
  for (l=1;l<=n_eqn;l++){
    phi(l,kp2+1) = phi(l,kp1+1);
    phi(l,kp1+1) = 0.0;
    p(l)       = 0.0;
  }
  for (int j=1;j<=k;j++){
    i     = kp1 - j;
    ip1   = i+1;
    temp2 = g(i+1);
    for (l=1;l<=n_eqn;l++){
      p(l)     = p(l) + temp2*phi(l,i+1);
      phi(l,i+1) = phi(l,i+1) + phi(l,ip1+1);
    }
  }
  if (nornd){
   p = y + p*h;}
   else{
    for (l=1;l<=n_eqn;l++){
      tau = h*p(l) - phi(l,16);
      p(l) = y(l) + tau;
      phi(l,17) = (p(l) - y(l)) - tau;
    }
  }
  xold = x;
  x = x + h;
  absh = fabs(h);
  yp = f(x,p);
  
  erkm2 = 0.0;
  erkm1 = 0.0;
  erk = 0.0;
  for (l=1;l<=n_eqn;l++){
    temp3 = 1.0/wt(l);
    temp4 = yp(l) - phi(l,1+1);
    if (km2> 0){
      erkm2 = erkm2 + ((phi(l,km1+1)+temp4)*temp3)*((phi(l,km1+1)+temp4)*temp3);
    }
    if (km2>=0){
      erkm1 = erkm1 + ((phi(l,k+1)+temp4)*temp3)*((phi(l,k+1)+temp4)*temp3);
    }
    erk = erk + (temp4*temp3)*(temp4*temp3);
  }
  
  if (km2> 0){
    erkm2 = absh*sig(km1+1)*gstr(km2+1)*sqrt(erkm2);
  }
  if (km2>=0){
    erkm1 = absh*sig(k+1)*gstr(km1+1)*sqrt(erkm1);
  }
  temp5 = absh*sqrt(erk);
  err = temp5*(g(k+1)-g(kp1+1));
  erk = temp5*sig(kp1+1)*gstr(k+1);
  knew = k;
  
  
  if (km2 >0){
    if (fmax(erkm1,erkm2)<=erk){
      knew=km1;
    }
  }
  if (km2==0){
    if (erkm1<=0.5*erk){
      knew=km1;
    }
  }
  
  success = (err<=epsilon);
  
  if (!success){
    phase1 = false; 
    x = xold;
    for (int i=1;i<=k;i++){
      temp1 = 1.0/beta(i+1);
      ip1 = i+1;
      for (l=1;l<=n_eqn;l++){
        phi(l,i+1)=temp1*(phi(l,i+1)-phi(l,ip1+1));
      }
    }
    
    if (k>=2){
      for (int i=2;i<=k;i++){
        psi_(i) = h-psi_(i+1);
      }
    }
    
    
    
    ifail = ifail+1;
    temp2 = 0.5;
    if (ifail>3) {
      if (p5eps < 0.25*erk){
        temp2 = sqrt(p5eps/erk);
      }
    }
    if (ifail>=3){
      knew = 1;
    }
    h = temp2*h;
    k = knew;
    if (fabs(h)<fouru*fabs(x)){
      crash = true;
      h = sign_(fouru*fabs(x), h);
      epsilon = epsilon*2.0;
        return y;                 
      }
  }  
  
  if (success){
    break;
  }
  
}

kold = k;
hold = h;

temp1 = h*g(kp1+1);
if (nornd){
  for (l=1;l<=n_eqn;l++){
    y(l) = p(l) + temp1*(yp(l) - phi(l,2));
  }
}
else{
  for (l=1;l<=n_eqn;l++){
    rhodouble = temp1*(yp(l) - phi(l,2)) - phi(l,17);
    y(l) = rhi + p(l);
    phi(l,16) = (y(l) - p(l)) - rhodouble;
  }
}
yp = f(x,y);


for (l=1;l<=n_eqn;l++){
  phi(l,kp1+1) = yp(l) - phi(l,2);
  phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
}
for (int i=1;i<=k;i++){
  for (l=1;l<=n_eqn;l++){
    phi(l,i+1) = phi(l,i+1) + phi(l,kp1+1);
  }
}


erkp1 = 0.0;
if ( (knew==km1) || (k==12) ){
  phase1 = false;
}

if (phase1){
  k = kp1;
  erk = erkp1;}
  else{
    if (knew==km1){
        
      k = km1;
      erk = erkm1;}
      else{
        if (kp1<=ns){
          for (l=1;l<=n_eqn;l++){
            erkp1 = erkp1 + (phi(l,kp2+1)/wt(l))*(phi(l,kp2+1)/wt(l));
          }
          erkp1 = absh*gstr(kp1+1)*sqrt(erkp1);
            
            
          if (k>1){
            if ( erkm1<=min(erk,erkp1)){
                    
              k=km1; erk=erkm1;}
              else{
                if ( (erkp1<erk) && (k!=12) ){
                        
                  k=kp1;
                  erk=erkp1;
                }
              }
            }
            else if (erkp1<0.5*erk){
                
              k = kp1;
              erk = erkp1;
            }
        } 
    } 
} 


if ( phase1 || (p5eps>=erk*two(k+2)) ){
  hnew = 2.0*h;}
  else{
    if (p5eps<erk){
      temp2 = k+1;
      r = pow(p5eps/erk,(1.0/temp2));
      hnew = absh*fmax(0.5, min(0.9,r));
      hnew = sign_(fmax(hnew, fouru*fabs(x)), h);}
      else{
        hnew = h;
      }
    }
    h = hnew;

    if (crash){
      State_    = DE_STATE.DE_BADACC;
      relerr    = epsilon*releps;       
      abserr    = epsilon*abseps;       
      y         = yy;                   
      t         = x;
      told      = t;
      OldPermit = true;
      return y;                       
    }

  nostep = nostep+1;  
  
  
  
  kle4 = kle4+1;
  if (kold>  4){
    kle4 = 0;
  }
  if (kle4>=50){
    stiff = true;
  }
}
/*
cout << "delsgn: " << delsgn << endl;
cout << "x: " << x << endl;
cout << "hi: " << hi << endl;
cout << "ki: " << ki << endl;
cout << "kold: " << kold << endl;
cout << "temp1: " << temp1 << endl;
cout << "term: " << term << endl;
cout << "psijm1: " << psijm1 << endl;
cout << "eta: " << eta << endl;
cout << "sum: " << sum << endl;
cout << "absh: " << absh << endl;
cout << "hold: " << hold << endl;
cout << "hnew: " << hnew << endl;
cout << "k: " << k << endl;
cout << "round: " << round << endl;
cout << "gamma: " << gamma << endl;
cout << "i: " << i << endl;
cout << "p5eps: " << p5eps << endl;
cout << "ifail: " << ifail << endl;
cout << "kp1: " << kp1 << endl;
cout << "kp2: " << kp2 << endl;
cout << "km1: " << km1 << endl;
cout << "km2: " << km2 << endl;
cout << "ns: " << ns << endl;
cout << "nsp1: " << nsp1 << endl;
cout << "realns: " << realns << endl;
cout << "im1: " << im1 << endl;
cout << "temp2: " << temp2 << endl;
cout << "temp3: " << temp3 << endl;
cout << "reali: " << reali << endl;
cout << "temp4: " << temp4 << endl;
cout << "nsm2: " << nsm2 << endl;
cout << "limit1: " << limit1 << endl;
cout << "temp5: " << temp5 << endl;
cout << "temp6: " << temp6 << endl;
cout << "limit2: " << limit2 << endl;
cout << "nsp2: " << nsp2 << endl;
cout << "ip1: " << ip1 << endl;
cout << "tau: " << tau << endl;
cout << "xold: " << xold << endl;
cout << "erkm2: " << erkm2 << endl;
cout << "erkm1: " << erkm1 << endl;
cout << "erk: " << erk << endl;
cout << "err: " << err << endl;
cout << "knew: " << knew << endl;
cout << "rhi: " << rhi << endl;
cout << "h: " << h << endl;
cout << "r: " << r << endl;
cout << "erkp1: " << erkp1 << endl;
cout << "rhodouble: " << rhodouble << endl;
cout << "told: " << told << endl;
cout << "epsilon: " << epsilon << endl;
cout << "del: " << del << endl;
cout << "absdel: " << absdel << endl;
cout << "tend: " << tend << endl;
cout << "nostep: " << nostep << endl;
cout << "kle4: " << kle4 << endl;
cout << "releps: " << releps << endl;
cout << "abseps: " << abseps << endl;
cout << "twou: " << twou << endl;
cout << "fouru: " << fouru << endl;

cout << "l: " << l << endl;
cout << "y: \n" << y << endl;*/

  cout<<69<<endl;            
return y;
}