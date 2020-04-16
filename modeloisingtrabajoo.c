#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define L 8
#define V L*L
#define N 1000
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
 int S[V];
 int xp[V],yp[V],xn[V],yyn[V];


FILE *fout;
FILE *fit;

double Random(void);

 void ini_ran(int SEMILLA);
void configuracionrandom(int *v);
void configuracioncgelada(int *w);
void medidas(double *ener, int *wv, double *mg);

void Histograma(double *data, double *Hist, int N_data, int N_Intervalos, double *d, double *m, double *M);
void calcularhistograma(double *X);
void vectorENMG( double *EE, double *MM, double energia, double magneti, int ii);
void generaconfiguracion(int flag);
void calcula_valores_medios(double *ee, double *mm, double *en, double *eren, double *mag, double *ermag, double *calor, double *xx,int i);
void inicivector( double *EE, double *MM,  int ii);
void generapn(int *xpo,int *ypo,int *xne,int *yne);
double error_medida_bloques(double *datos, int n_datos, int n_intervalos);

double media(double *serie, int numero);
double varianza(double *serie, int numero);
double cuadradosmedia(double *serie, int numero);


void metropolis(int *s,int *ntot, double beta);
int main()
{
        ini_ran(123456789);
 double beta_inicial, beta_final, delta_beta,beta, ener, mg;
 double E[1000],M[1000];
 double energia,error_energia,magnetizacion,error_magnetizacion,calor_esp,susceptibilidad;
 int N_ter, N_med, N_Met, N_pasos, N_betas,sentido,nacepta,N_M,N_m,i,ik;

 beta_inicial=0.40;
 beta_final=0.47;
 delta_beta=0.001;
 N_ter=200;
 N_med=1000;
 N_Met=5;
 generaconfiguracion(1);


 N_pasos=(beta_final-beta_inicial)/delta_beta;

 beta=beta_inicial;
 generapn(xp,yp,xn,yyn);



  fout=fopen("errires.txt","wt");
 for(sentido=0;sentido<1;sentido++)
 {
     for(N_betas=0;N_betas<N_pasos;N_betas++)
     {
         for(N_Met=0;N_Met<N_ter;N_Met++)
         {


              metropolis(S,&nacepta, beta);
            // medidas(&ener, S,&mg);
            //fprintf(fout," %d %f \n",N_Met,ener);

        }
     //exit(1);

        inicivector(E,M,N_med);

        for(N_M=0;N_M<N_med;N_M++)
        {
            for(N_m=0;N_m<N_Met;N_m++)
            {
                metropolis(S,&nacepta, beta);
            }
            medidas(&ener, S,&mg);
            //printf("%lf %lf\n",ener,mg);getchar();


             E[N_M]=ener;
             M[N_M]=mg;
                 //printf(" %f %f\n",E[N_M],M[N_M]);

        }

         //for(ik=0;ik<N_med;ik++)
           {
               //printf(" %f\n",E[ik]);
           }





                     //calcularhistograma(E);


    calcula_valores_medios(E, M, &energia, &error_energia,&magnetizacion,&error_magnetizacion, &calor_esp, &susceptibilidad,N_med);
      //printf(" %f \n",magnetizacion);



   //fprintf(fout," %f %f %f %f %f \n",beta,calor_esp,susceptibilidad, error_energia, error_magnetizacion);
     fprintf(fout," %f %f  \n",beta,calor_esp);



          beta+=delta_beta;

      }



            delta_beta=-delta_beta;
     }


 fclose(fout);

}







void generaconfiguracion(int flag)
{
    switch(flag)
    {
        case 1:
         configuracionrandom(S);
         break;
        case 0:
            configuracioncgelada(S);
            break;
        case -1:
            break;


    }
}
void vectorENMG( double *EE, double *MM, double energia, double magneti, int ii)
{

   EE[ii]=energia;
   MM[ii]=magneti;


}
void inicivector( double *EE, double *MM,  int ii)
{
    for(int i=0;i<ii;i++)
    {
        EE[i]=0;
        MM[i]=0;
    }
}
void calcularhistograma(double *X)
{
    int j,i;
        fout=fopen("hista.dat","wt");

    double H[10000],delta,m,n;

    Histograma(X,H,N,100,&delta,&n,&m);
     for(j=0;j<100;j++)
      {

          fprintf(fout," %d %f %f\n",j,j*delta+n,H[j]);


      }

    fclose(fout);

}
void Histograma(double *data, double *Hist, int N_data, int N_Intervalos, double *d, double *m, double *M)
 {
     int i,Indice;
     double del,min,max,Norm;
     for(i=0;i<N_Intervalos;i++)
     {
         Hist[i]=0;
     }
     min=10000000;
     max=-10000000;

     for(i=0;i<N_data;i++)
     {
         if(data[i]<min)min=data[i];
         if(data[i]>max)max=data[i];

     }

     del=(max-min)/N_Intervalos;
     if(del==0)
     {
         printf("error: no se pueden calcular los intervalos");
         exit(1);
     }


     for(i=0;i<N_data;i++)
     {
         Indice=(data[i]-min)/del;
         Hist[Indice]++;
     }
     *d=del;
     *m=min;
     *M=max;

     Norm=1.0/(N_data*del);
     for(i=0;i<N_Intervalos;i++)
     {
         Hist[i]*=Norm;
     }
 }
void configuracionrandom(int *v)
{
    int i;
    for(i=0;i<V;i++)
    {


            if(Random()<0.5)
                v[i]=1;
                else
                v[i]=-1;

    }


}
void configuracioncgelada(int *w)
{
    int i;
    for(i=0;i<V;i++)
    {


                w[i]=1;
                 w[i+1]=-1;

    //printf(" %d \n",i);
    i=1+i;
    }

}
double Random(void)
{
    double r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}
void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR*SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

void medidas(double *ener, int *wv, double *mg)
{
    //int xp[L],yp[L];
    int i,j,n;
   /* for(i=0;i<(L-1);i++)
    {
        xp[i]=1;
        xp[L-1]=-(L-1);
    }
    for(i=0;i<(L-1);i++)
    {
        yp[i]=L;
        yp[L-1]=-L*(L-1);
    }*/
     int sum,summ;
    sum=0;
    n=0;
    for(j=0;j<L;j++)
    {
        for(i=0;i<L;i++)
        {

            sum+=-wv[n]*(wv[n+xp[i]]+wv[n+yp[j]]);
                // printf(" %d %d %d %d \n",sum,wv[n],wv[n+xp[i]],wv[n+yp[j]]);

              n++;



        }
    }
     summ=0;
    for(i=0;i<V;i++)
    {
        summ+=wv[i];
    }
    *mg=fabs(summ/((double)V));


    *ener=sum/(2*((double)V));
}






void metropolis(int *s,int *ntot, double beta)
{
    int x,y,i,ind,n,in;
    double prob[5];
   prob[0]=exp(-beta*(-8));
   prob[1]=exp(-beta*(-4));
   prob[2]=exp(-beta*0);
   prob[3]=exp(-beta*4);
   prob[4]=exp(-beta*8);

   n=0;

   for(y=0;y<L;y++)
   {
      for(x=0;x<L;x++)
   {
       ind=s[n]*(s[n+xp[x]]+s[n+yp[y]]+s[n+xn[x]]+s[n+yyn[y]])/2+2;

             // for(i=0;i<L;i++)
              //printf(" %d %d %d %d %d %d\n",ind,s[n],s[n+xp[x]],s[n+yp[y]],s[n+xn[x]],s[n+yyn[y]]);

       if(Random()<prob[ind])
       {
             s[n]=-s[n];
       }
        n++;
   }

   }


    *ntot=n;


}
void generapn(int *xpo,int *ypo,int *xne,int *yne)
{
    int i;
     for(i=0;i<(L-1);i++)
    {
        xpo[i]=1;
        xpo[L-1]=-(L-1);
    }
    for(i=0;i<(L-1);i++)
    {
        ypo[i]=L;
        ypo[L-1]=-L*(L-1);
    }
     for(i=1;i<L;i++)
    {
        xne[i]=-1;
    }
    xne[0]=(L-1);

    for(i=1;i<L;i++)
    {
        yne[i]=-L;
    }
    yne[0]=L*(L-1);


}
void calcula_valores_medios(double *ee, double *mm, double *en, double *eren, double *mag, double *ermag, double *calor, double *xx,int i)
{
    double sum,sumv,sumva,sumvaa,eE,mM,enn,emm;
    int j;
       sum=0;
       sumv=0;
       sumva=0;
       sumvaa=0;
       *eren=error_medida_bloques(ee, i, 100);
       *ermag=error_medida_bloques(mm, i, 100);
     for(j=0;j<i;j++)
     {
       sum+=ee[j];
      }
      *en=sum/((double)i);
       enn=sum/((double)i);

     for(j=0;j<i;j++)
     {
       sumv+=mm[j];
      }
      *mag=sumv/((double)i);
      emm=sumv/((double)i);

      for(j=0;j<i;j++)
      {
          sumva+=ee[j]*ee[j];

      }
    eE=sumva/((double)(i*i));
    for(j=0;j<i;j++)
      {
          sumvaa+=mm[j]*mm[j];

      }
    mM=sumvaa/((double)(i*i));
      enn=2*V*varianza(ee, i);
    *calor=enn;
//    *xx=V*((emm*emm)-((sumv/i)*(sumv/i)));
}
double error_medida_bloques(double *datos, int n_datos, int n_intervalos)
{
    int interv_size=n_datos/n_intervalos,i,j;
    double datosauxiliar[interv_size],mediasint[n_intervalos];
    double varianza_int;

    for(j=0;j<n_intervalos;j++)
    {
        for(i=0;i<interv_size;i++)
            datosauxiliar[i]=datos[j*interv_size+i];


            mediasint[j]=media(datosauxiliar, interv_size);

    }
    varianza_int=varianza(mediasint,n_intervalos);



    return sqrt(varianza_int/interv_size);
}
double media(double *serie, int numero)
{
    int i,j;
    double x=0;
    for(i=0;i<numero;i++)
    {
        x+=serie[i];
    }
    return x/numero;


}
double cuadradosmedia(double *serie, int numero)
{
    int j;
    double y;
    y=0;
    for(j=0;j<numero;j++)
    {


        y+=(serie[j]*serie[j]);
         //printf(" %f\n",y);
    }


    return y/numero;

}
double varianza(double *serie, int numero)
{
    int i,j;
    double z;
    z=0;
    double mediacuadrados,Media;
    mediacuadrados=cuadradosmedia(serie,numero);


    Media=media(serie,numero);
    z=(mediacuadrados-(Media*Media));
    return z;

}
