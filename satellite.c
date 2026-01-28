#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1e5
#define G 6.67e-11
#define M 6e24

int main(){

  int i=1,choice;
  double v,dv,r=4e5/*5*/,dr,t=0,dt=0.01,dm,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,l=1,m=6000,mmin,E,dE,E0,E1;
  FILE* fp;

  if((fp=fopen("sat.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
  
  v = sqrt(G*M/r);
  //v=0;
  E0 = G*M/r + 0.5*m*v*v;
  E = 1;
  mmin = 5000;
  printf("1\t2\t3\n");
  scanf("%d",&choice);

  if(choice == 1){

  do{

    if(E<=0.000000001/*r<=100*/){
      i=0;
    }

    fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",t,r,v,E);
     
    a1 = 4*(l*l*l*l*(r)/(m*m))*dt;
    b1 = -2*((l*l*(r))/m)*dt;
    
    a2 = 4*(l*l*l*l*(r+0.5*b1)/(m*m))*dt;
    b2 = -2*((l*l*(r+0.5*b1))/m)*dt;
    
    a3 = 4*(l*l*l*l*(r+0.5*b2)/(m*m))*dt;
    b3 = -2*((l*l*(r+0.5*b2))/m)*dt;
    
    a4 = 4*(l*l*l*l*(r+b3)/(m*m))*dt;
    b4 = -2*((l*l*(r+b3))/m)*dt;

    dv = (a1+2*a2+2*a3+a4)/6;
    dr = (b1+2*b2+2*b3+b4)/6;

    v+=dv;
    r+=dr;
    t+=dt;
    E += -(l*l*G*M/(m*r))/E0;
    
    
    
  } while(i!=0); 
  printf("\n%lf\t%lf\t%lf\t%lf\n\n",t,r,v,E);
  }
  
  i = 1;             // Reimpostazione parametro di break

  E = E0;
  
  if(choice == 2){

    do{

      if(E<=0.01){
	i=0;
      }
      
      if(m > mmin){
	
	fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\tZ\n",t,r,v,m,E);
	
	c1 = -(dt*G*M*(m)/((v)*(r)*(r))) - (l*l*dt);
	a1 = -dt*G*M/((r)*(r)) - ((v)/(m))*(l*l*dt+c1);
	b1 = -2*((l*l*(r))/(m))*dt - 2*((r)/(m))*c1;
	
	c2 = -(dt*G*M*(m+0.5*c1)/((v+0.5*a1)*(r+0.5*b1)*(r+0.5*b1))) - (l*l*dt);
	a2 = -dt*G*M/((r+0.5*b1)*(r+0.5*b1)) - ((v+0.5*a1)/(m+0.5*c1))*(l*l*dt+c2);
	b2 = -2*((l*l*(r+0.5*b1))/(m+0.5*c1))*dt - 2*((r+0.5*b1)/(m+0.5*c1))*c2;
	
	c3 = -(dt*G*M*(m+0.5*c2)/((v+0.5*a2)*(r+0.5*b2)*(r+0.5*b2))) - (l*l*dt);
	a3 = -dt*G*M/((r+0.5*b2)*(r+0.5*b2)) - ((v+0.5*a2)/(m+0.5*c2))*(l*l*dt+c3);
	b3 = -2*((l*l*(r+0.5*b2))/(m+0.5*c2))*dt - 2*((r+0.5*b2)/(m+0.5*c2))*c3;
	
	c4 = -(dt*G*M*(m+c3)/((v+a3)*(r+b3)*(r+b3))) - (l*l*dt);
	a4 = -dt*G*M/((r+b3)*(r+b3)) - ((v+a3)/(m+c3))*(l*l*dt+c4);
	b4 = -2*((l*l*(r+b3))/(m+c3))*dt - 2*((r+b3)/(m+c3))*c4;
	
	dv = (a1+2*a2+2*a3+a4)/6;
	dr = (b1+2*b2+2*b3+b4)/6;
	dm = (c1+2*c2+2*c3+c4)/6;
	
	v+=dv;
	r+=dr;
	m+=dm;
	t+=dt;
	
	E += -v*v*(l*l*dt+dm);
	
	
	
      } else {

	fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",t,r,v,m,E);
	
	c1 = 0;
	c2 = 0;
	c3 = 0;
	c4 = 0;
	dm = 0;
	
	a1 = 4*(l*l*l*l*(r)/(m*m))*dt;
	b1 = -2*((l*l*(r))/m)*dt;
	
	a2 = 4*(l*l*l*l*(r+0.5*b1)/(m*m))*dt;
	b2 = -2*((l*l*(r+0.5*b1))/m)*dt;
	
	a3 = 4*(l*l*l*l*(r+0.5*b2)/(m*m))*dt;
	b3 = -2*((l*l*(r+0.5*b2))/m)*dt;
	
	a4 = 4*(l*l*l*l*(r+b3)/(m*m))*dt;
	b4 = -2*((l*l*(r+b3))/m)*dt;
	
	dv = (a1+2*a2+2*a3+a4)/6;
	dr = (b1+2*b2+2*b3+b4)/6;
	
	v+=dv;
	r+=dr;
	t+=dt;
	
	E += -(l*l*G*M/(m*r));
	
	
	
      }
    } while(i!=0); 
    printf("\n%lf\t%lf\t%lf\t%lf\n\n",t,r,v,E);
  }



   if(choice == 3){

    do{

      if(E<=0.01){
	i=0;
      }

       if(m > mmin){
	
	fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\tZ\n",t,r,v,m,E);

	m += -l*l*dt;
	v -= -G*M*dt/(r*r*10000);
	t += dt;
	  
       } else {

	 fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\tZ\n",t,r,v,m,E);
	 
	 a1 = 4*(l*l*l*l*(r)/(m*m))*dt;
	 b1 = -2*((l*l*(r))/m)*dt;
	 
	 a2 = 4*(l*l*l*l*(r+0.5*b1)/(m*m))*dt;
	 b2 = -2*((l*l*(r+0.5*b1))/m)*dt;
	 
	 a3 = 4*(l*l*l*l*(r+0.5*b2)/(m*m))*dt;
	 b3 = -2*((l*l*(r+0.5*b2))/m)*dt;
	 
	 a4 = 4*(l*l*l*l*(r+b3)/(m*m))*dt;
	 b4 = -2*((l*l*(r+b3))/m)*dt;
	 
	 dv = (a1+2*a2+2*a3+a4)/6;
	 dr = (b1+2*b2+2*b3+b4)/6;
	 
	 v+=dv;
	 r+=dr;
	 t+=dt;
	 
	 E += -(l*l*G*M/(m*r));
	 
       }
       
    } while(i!=0);
    
    
   }
    printf("\n%lf\t%lf\t%lf\t%lf\n\n",t,r,v,E);
}



///////////////////////////////
//////  SCRIVERE dr = 0 ///////
///////////////////////////////
