#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define N 10
//#define DEBUG
//#define DEBUG3
//#define DEBUG4
#define DEBUG5
void configuracion_puntos_plano(double *x,double *y,int flag);

double Norma(double *p);

void configuracion_inicial(int *K,int flag);

double Energia(double *x,double *y,int *K);

void cambio_tentativo_A(int *K,int *K_aux);

unsigned int irr[256];

unsigned int ir1;

unsigned char ind_ran,ig1,ig2,ig3;

double Random(void);

void ini_ran(int SEMILLA);

void copia(double *p,double *q);

void inicia_vector(double *p,int n,int flag);

FILE *fout;
FILE *aceptan;

int main()
{
    int K[N+1];
    int K_aux[N+1];
    ini_ran(1234567);
    int i,ITER,INER,iter,iner,tot,acept;
    double x[N+1],y[N+1];
    double beta,epsilon,delta_beta,delta_epsilon,aceptancia;
    double Resto,Resto_n;
    double aux1,aux2;
    double Ener;
    //Inicializamos los vectores
    x[0]=y[0]=K[0]=K_aux[0];

    fout=fopen("evol_viajante.txt","wt");

    aceptan=fopen("aceptancia_viajante.txt","wt");

    beta=1.0*sqrt((double)N);
    epsilon=0.1/sqrt((double)N);

    configuracion_puntos_plano(x,y,0);
    configuracion_inicial(K,1);


    ITER=2;
    INER=2;
    delta_epsilon=-epsilon/(ITER+1);
    delta_beta=beta/ITER;

    for(iter=0;iter<ITER;iter++)
    {
                                    tot=acept=0;
                                    beta+=delta_beta;
                                    epsilon+=delta_epsilon;

                    for(iner=0;iner<INER;iner++)
                        {
                            Resto=Energia(x,y,K);

                           // cambio_tentativo_A(K,K_aux);
                            cambio_tentativo_B(K,K_aux);
                            Resto_n=Energia(x,y,K_aux);
                            Ener=beta*(Resto_n-Resto);
                            #ifdef DEBUG
                            {
                                for(i=0;i<N+1;i++)
                                    {
                                    printf("x[%d]:%lf, y[%d]:%lf,K[%d]:%d,K_aux[%d]:%d,Ener:%lf\n",i,x[i],i,y[i],i,K[i],i,K_aux[i],Ener);
                                    }
                                printf("\n\n");
                            }
                            #endif // DEBUG
                           // printf("Ener:%lf\n",Ener);
                            tot++;
                           if(exp(-Ener)>Random())
                                {

                                    copia(K_aux,K);
                                    Resto=Resto_n;
                                    acept++;
                                }

                            fprintf(fout,"%lf %d %d %lf %lf\n",beta,tot,acept,log(Resto),log(Resto_n));
                        }
                       // printf("\n\n\nhola!\n\n\n");
            aceptancia=(double)acept/((double)tot);
            fprintf(aceptan,"%lf %lf\n",aceptancia,beta);
    }


    #ifdef DEBUG
    {
       for(i=0;i<N+1;i++)
            {
                printf("x[%d]:%lf, y[%d]:%lf,K[%d]:%d,K_aux[%d]:%d\n",i,x[i],i,y[i],i,K[i],i,K_aux[i]);
            }

           printf("x[0]=%lf,y[0]=%lf,x[N+1]=%lf,y[N+1]=%lf,Resto:%lf\n\n",x[0],y[0],x[N],y[N],Resto);

    }
    #endif // DEBUG

                    fclose(fout);
                    fclose(aceptan);

    return 0;
}

double Energia(double *x,double *y,int *K)
{
    int i;

    double sum,aux;

    sum=0;

    for(i=0;i<N;i++)
        {
            sum+=(x[K[i]]-x[K[i+1]])*(x[K[i]]-x[K[i+1]])+(y[K[i]]-y[K[i+1]])*(y[K[i]]-y[K[i+1]]);

            #ifdef DEBUG4
            aux=(x[K[i]]-x[K[i+1]])*(x[K[i]]-x[K[i+1]])+(y[K[i]]-y[K[i+1]])*(y[K[i]]-y[K[i+1]]);
            printf("suma parcial:%lf\n",aux);
            #endif // DEBUG4
        }

    sum=sqrt(sum);

    return sum;
}

void copia(double *p,double *q)
{
    int i;

    for(i=0;i<N;i++)
        q[i]=p[i];
        q[N]=p[0];
}

void configuracion_puntos_plano(double *x,double *y,int flag)
{
        int i;
        double aux1,aux2;

        if(flag==0) //puntos aleatorios
        {
           for(i=0;i<N;i++)

                {
                  x[i]=(double)Random();
                  y[i]=(double)Random();
                }

                x[N]=x[0];
                y[N]=y[0];

        }


}

void configuracion_inicial(int *K,int flag)
{
    int i,j,aux1,aux2;

    if(flag==0) //Configuracion random
    {
        for(i=0;i<N;i++)
        {
            K[i]=i;
        }
        for(i=0;i<N;i++)
        {
           j=(int)(N*Random());
           aux1=K[i];
           aux2=K[j];
           K[i]=aux2;
           K[j]=aux1;
        }
        K[N]=K[0];

    }

    if(flag==1) //Configuracion permutacion zero
    {
        for(i=0;i<N;i++)
        {
            K[i]=i;
        }
        K[N]=K[0];
    }
    /*
    if(flag==2) //Configuracion del tonto
    {
        for(i=0;i<N;i++)
        {

        }

    }
    */
}
//void beta_inicial()
void cambio_tentativo_A(int *K,int *K_aux)
{
    int ini,fin,i,j,aux1,aux2;
    ini=(int)(N*Random());
    fin=(int)(N*Random());
    #ifdef DEBUG3
    int auxextra;
    auxextra=ini;
    #endif // DEBUG3

    j=0;


    if(ini==fin)
        while(ini==fin)
            ini=(int)(N*Random());

    if(ini>fin)
    {
        aux1=ini;
        aux2=fin;
        ini=aux2;
        fin=aux1;
        #ifdef DEBUG3
        auxextra=ini;
        #endif // DEBUG3
    }

   for(i=0;i<N;i++)
        K_aux[i]=K[i];
    K_aux[N]=K[0];

   for(i=ini;i<=fin;i++)
   {
       K_aux[i]=K[fin-j];
       j++;
   }
   K_aux[N]=K_aux[0];

    #ifdef DEBUG3
    for(j=0;j<N;j++)
            printf("K[%d]:%d\n",j,K[j]);
            printf("cambio entre la posicion %d y la posicion %d\n",auxextra,fin);
    #endif // DEBUG3



}
void cambio_tentativo_B(int *K,int *K_aux)
{
    int a,b,delta,aux1,aux2,i,j,dif;
    int intervalo1,intervalo2;
    a=(int)(N*Random()); // lugar donde insertamos el cambio
    b=(int)(N*Random());	// lugar donde cojemos el cambio
    delta=(int)(N/2*Random()); // intervalo que cojemos de numeros con los que haremos el cambio. Está ajustado entre N/2 para que no nos de problemas

    if(a==b)
        while(a==b)
            a=(int)(N*Random());

    if(a>b)
    {
        aux1=a;
        aux2=b;
        a=aux2;
        b=aux1;
    }
  
	while(b+delta<N)
		delta=(int)(N/2*Random()); // este bucle sirve para que el delta no se nos salga del intervalo

    for(i=0;i<N;i++)
    K_aux[i]=K[i];
    K_aux[N]=K[0];
    
   for(i=a;i<=a+delta;i++)
   {
       K_aux[i]=K[b+i];

   }
	for(i=a+delta+1;i<=N-a;i++)
	{
		
		K_aux[i]=K[i-delta];
	}


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

//printf("r=%f\n",r);

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

INI=(INI*FACTOR+SUM);

irr[i]=INI;

}

ind_ran=ig1=ig2=ig3=0;

}


