#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//#define DEBUG
#define Dim_Max 2
#define N_par_Max 6000
#define NLinks_Max (N_par_Max*200)
#define fran ((double)rand()/((double)RAND_MAX+1.0))

void lee_param(int ac,char **av);
void lee_coords(void);
void genera_coords(void);
void lee_graph(void);
void genera_graph(void);
void Calcula_Potential_FR(void);
void escribe_coords(void);
void compute_degree(void);
void Evoluciona(void);
void Acumula_Links(void);

double usada[N_par_Max][N_par_Max];
double F[N_par_Max][Dim_Max];
double epsilon;

struct Prm
{
            int Dim;
            int N_par;
            int niter;
            int flag_graph;
            int flag_coords;
            char graph_file[256];
            char coords_file[256];
}
Param;

struct grph
{
            int N_par;
            int N_links;
            double wmean;
            double grado[N_par_Max];
            double size[N_par_Max];
            double r[N_par_Max][Dim_Max];
            int origen[NLinks_Max];
            int destino[NLinks_Max];
            double peso[NLinks_Max];
}
graph;
struct P_FR
{
            double area;
            double maxdelta;
            double coolexp;
            double temperature;
            double repulserad;
            double frk;
            double frk2;
}
Potential_FR;

int Lado_Max;

int main(int argc,char **argv)
{
                int i,segundos;

                int mesfr;

                FILE *fevol;

                epsilon=0.000001;
                 //para evitar la singularidad en el origen

                lee_param(argc,argv);

                if(Param.flag_graph==0)
                    genera_graph();

                else
                    lee_graph();

                compute_degree();

                Lado_Max=sqrt((double)graph.N_par);

                if(Param.flag_coords==0)
                    genera_coords();

                else
                    lee_coords();

                Calcula_Potential_FR();

                mesfr=10;

                fevol=fopen("Evol.dat","wt");

                segundos=time(0);

                for (i=Param.niter;i>0;i--)
                {
                    Potential_FR.temperature=Potential_FR.maxdelta*pow(i/(double)Param.niter,Potential_FR.coolexp);
                    Evoluciona();
                    fprintf(fevol," %lf %lf %lf %lf\n",graph.r[0][0],graph.r[0][1],graph.r[1][0],graph.r[1][1]);

                    if(i%mesfr==0)
                    {
                        printf("(t=%d)",(int)(time(0)-segundos));
                        printf(" %d %lf %lf %lf %lf \n",i,graph.r[0][0],graph.r[0][1],graph.r[1][0],graph.r[1][1]);
                        //getchar();
                    }
                }

                printf("(Total Time= %d\n",(int)(time(0)-segundos));

                fclose(fevol);

                escribe_coords();

                return 1;
}

void Evoluciona(void)
{
    double distance,d2;

    int i,j,il,k;

    double factor,del[Dim_Max],ded;

    //Inicializamos el vector doble
    for(i=0;i<graph.N_par;i++)
        for(k=0;k<Param.Dim;k++)
            F[i][k]=0;

    for(i=0;i<graph.N_par;i++)

        {
            for(j=i+1;j<graph.N_par;j++)

                {
                    for(k=0,d2=0;k<Param.Dim;k++)
                        {
                            del[k]=graph.r[i][k]-graph.r[j][k];
                            d2+=del[k]*del[k];
                        }

                d2+=epsilon;

                distance=sqrt(d2);

                factor=Potential_FR.frk2*(1.0/d2-distance/Potential_FR.repulserad);

                for(k=0;k<Param.Dim;k++)
                    {
                        F[i][k]+=del[k]*factor;
                        F[j][k]-=del[k]*factor;
                    }
                }
        }

    for(il=0;il<graph.N_links;il++) //recorre subconjunto de links
    {
        i=graph.origen[il];

        j=graph.destino[il];

        for(k=0,d2=0;k<Param.Dim;k++)
            {
                del[k]=graph.r[i][k]-graph.r[j][k];
                d2+=del[k]*del[k];
            }

        distance=sqrt(epsilon+d2);


        factor=distance/Potential_FR.frk*graph.peso[il];

        for(k=0;k<Param.Dim;k++)
            F[i][k]-=del[k]*factor;

        F[j][k]+=del[k]*factor;
    } //fin bucle atraccion
    for(i=0;i<graph.N_par;i++) //recorre subconjunto de nodos
    {
        for(k=0,ded=0;k<Param.Dim;k++)
            ded+=F[i][k]*F[i][k];

        ded=sqrt(ded+epsilon);

        if(ded>Potential_FR.temperature) //Una especie de Simmulated Annealing
            {
                ded=Potential_FR.temperature/ded;

                for(k=0;k<Param.Dim;k++)
                    F[i][k]*=ded;
            }

        for(k=0;k<Param.Dim;k++)
            graph.r[i][k]+=F[i][k]; //Consumado el cambio
    } //fin bucle acumula y cambia
} //fin de funcion

void lee_param(int ac,char **av)
{
char dummy[256],name[256];
FILE *Dummy_File,*Input_File;
/*
Leemos el fichero que se indica en linea de comandos (parameters.dat default);
En ese fichero leemos parametros para la simulacion
*/
Param.Dim=2;

switch(ac)
    {
        case 1:
        sprintf(name,"parameters.dat");
        printf("Opening parameters.dat for input\n");
        break;

        case 2:
        sscanf(av[1],"%s",name);
        break;

        default:
        printf("uso: fr_new2 [parameters_file:Def:paramemeters.dat] \n");
        exit(1);
    }

    if((Input_File=fopen(name,"rt"))==NULL)
    {
        printf("Error opening file %s for reading\n",name);
        exit(1);
    }

    if(fscanf(Input_File,"%d%s\n",&Param.N_par,dummy)!=2)
    {
        printf("Error reading first data N_par=%d\n",Param.N_par);
        exit(2);
    }

    if(fscanf(Input_File,"%d%s\n",&Param.niter,dummy)!=2)
    {
        printf("Error reading data niter=%d\n",Param.niter);
        exit(3);
    }

    if(fscanf(Input_File,"%d%s\n",&Param.flag_graph,Param.graph_file)!=2)
    {
        printf("Error reading parameters for graph_file\n");
        exit(4);
    }

    else
    {
        if(Param.flag_graph==1)
        {
            if((Dummy_File=fopen(Param.graph_file,"rt"))==NULL)
            {
                printf("Error opening Links File=%s\n",Param.graph_file);
                exit(5);
            }
            else
                fclose(Dummy_File);
        }
    }

    if(fscanf(Input_File,"%d%s\n",&Param.flag_coords,Param.coords_file)!=2)
    {
        printf("Error reading parameters for coords_file\n");
        exit(6);
    }

    else
    {
        if(Param.flag_coords==1)
        {
            if((Dummy_File=fopen(Param.coords_file,"rt"))==NULL)
            {
                printf("Error opening Coords File=%s\n",Param.coords_file);
                exit(7);
            }

            else
                fclose(Dummy_File);
        }
    }

    fclose(Input_File);
    #ifdef DEBUG
    printf("Leidos Parametros en lee_param:\n");
    printf("N_par=%d\n",Param.N_par);
    printf("flag_graph=%d\n",Param.flag_graph);
    printf("flag_coords=%d\n",Param.flag_coords);
    #endif
}

void lee_coords(void)
{
    FILE *fin;

    int site;

    if((fin=fopen(Param.coords_file,"rt"))==NULL)
    {
        printf("Error en apertura de fichero de coordenadas: %s\n",Param.coords_file);
        exit(30);
    }

    site=0;

    while(fscanf(fin,"%lf %lf\n",&graph.r[site][0],&graph.r[site][1])!=EOF)site++;
        if(site!=Param.N_par)
            {
                printf("Number of particles modified reading coords_file:\n");
                printf("N in Input_file=%d,N in coords_file=%d\n",Param.N_par,site);
                printf("All OK\n");
            }

    fclose(fin);

    Param.N_par=site;

    #ifdef DEBUG
    printf("Coordenada r[0][0]=%lf, r[%d][1]=%lf\n",
    graph.r[0][0],graph.N_par-1,graph.r[graph.N_par-1][0]);
    #endif
}

void genera_coords(void)
{
    int i,j;
    for(i=0;i<graph.N_par;i++)
        for(j=0;j<Param.Dim;j++)
            graph.r[i][j]=Lado_Max*(fran-0.5);
    for(j=0;j<Param.Dim;j++)
        graph.r[0][j]=0;
}

void lee_graph(void)
{
    //Formato del archivo: n_ini n_fin link_ini_fin
    //n_ini debe ser un indice de 0 a (n_par-1)
    int ilink,np,A,B;

    FILE *fin;
    graph.wmean=0;
    ilink=np=0;

    fin=fopen(Param.graph_file,"rt");

    while( fscanf(fin,"%d %d %lf\n",&A,&B,&graph.peso[ilink])!= EOF)
    {
        if(A!=B)//NO queremos loops
        {
            graph.origen[ilink]=A;
            graph.destino[ilink]=B;
            graph.wmean+=graph.peso[ilink];
            #ifdef DEBUG
            printf("site_i=%d,site_j=%d,link[%d]=%lf\n",A,B,ilink,graph.peso[ilink]);
            #endif
            ilink++;

            if(A>np)np=A;

            if(B>np)np=B;
        }
    }
    if((np+1)!=Param.N_par)
        printf("Sobreescrito el numero de particulas: input=%d, graph_file=%d\n",Param.N_par,np+1);

    graph.N_par=np+1; //se sobrescribe lo leido en el fichero de parametros.

    graph.N_links=ilink;

    fclose(fin);

    if(graph.N_links>0)
        graph.wmean/=ilink;

    if(graph.N_par<=0)
            {
                printf("Error: Grafo leido tiene 0 nodos\n");
                exit(10);
            }
    #ifdef DEBUG
    printf("Grafo:\n");
    printf("N_par=%d\n",graph.N_par);
    printf("N_links=%d\n",graph.N_links);
    printf("wmean=%lf\n",graph.wmean);
    #endif
}
void genera_graph(void)
{
    int N_Links_Medio,ilink,i,j;
    int A,B;
    double Peso_Maximo_Link;
    graph.N_par=Param.N_par;
    ilink=0;
    for(i=0;i<graph.N_par;i++)
    for(j=0;j<graph.N_par;j++)
    usada[i][j]=0;
    N_Links_Medio=(sqrt(sqrt(graph.N_par)));

    if(N_Links_Medio<2)N_Links_Medio=2;

    if(Param.N_par>=1000)N_Links_Medio=1;

    Peso_Maximo_Link=20.0;
    for(i=0;i<graph.N_par*N_Links_Medio;i++) //random
    {
        do
            {
                A=graph.N_par*fran;

                while((B=graph.N_par*fran)==A);

            }

        while(usada[A][B]==1);

        graph.origen[ilink]=A;

        graph.destino[ilink]=B;

        usada[A][B]=1;

        graph.peso[ilink]=Peso_Maximo_Link*fran;

        graph.wmean+=graph.peso[ilink];

        ilink++;
    }
    graph.N_links=ilink;
    if(graph.N_links>0)
    graph.wmean/=ilink;
    if(graph.N_par<=0)
    {
    printf("Error: Grafo construido tiene 0 nodos\n");
    exit(10);
    }
}

void Calcula_Potential_FR(void)
{
    Potential_FR.area=(double)graph.N_par*graph.N_par;
    Potential_FR.maxdelta=sqrt(Potential_FR.area);
    Potential_FR.coolexp=1.5;
    Potential_FR.repulserad=graph.N_par*Potential_FR.area;
    Potential_FR.frk=sqrt((double)graph.N_par)*pow(graph.wmean,1/3.0);
    Potential_FR.frk2=Potential_FR.frk*Potential_FR.frk;
    #ifdef DEBUG
    printf("En Calcula_Potential_FR:\narea=%lf, maxdelta=%f,repulserad=%lf, frk=%lf\n",
    Potential_FR.area,Potential_FR.maxdelta,Potential_FR.repulserad,Potential_FR.frk);
    #endif
}

void escribe_coords(void)
{
            FILE *fout;
            int i,j,il,k;
            fout=fopen("posicion.dat","wt");
            for(i=0;i<graph.N_par;i++)
                {
                    for(k=0;k<(Param.Dim-1);k++)
                        fprintf(fout,"%lf ",graph.r[i][k]);
                    fprintf(fout,"%lf\n",graph.r[i][Param.Dim-1]);//para evitar el blanco al final de linea
                }
            fclose(fout);

            fout=fopen("posicion_r.dat","wt");

            for(i=0;i<graph.N_par;i++)
            {
                for(k=0;k<Param.Dim;k++)
                    fprintf(fout,"%lf ",graph.r[i][k]);
                fprintf(fout," %lf\n",sqrt(1+graph.grado[i]));
            }
            fclose(fout);

            fout=fopen("posicion.plt","wt");

            if(Param.Dim==2)
                fprintf(fout,"plot \"posicion_r.dat\" u 1:2:%d with points pt 6 ps variable\n",
                Param.Dim+1);

            else
                fprintf(fout,"splot \"posicion_r.dat\" u 1:2:3:%d with points pt 6 ps variable\n",
                Param.Dim+1);

            fprintf(fout,"unset arrow\n");

            for(il=0;il<graph.N_links;il++)
            {
                i=graph.origen[il];
                j=graph.destino[il];
                fprintf(fout,"set arrow from ");
                for(k=0;k<Param.Dim-1;k++)
                    fprintf(fout,"%014.8f,\t",graph.r[i][k]);

                fprintf(fout,"%014.8f ",graph.r[i][Param.Dim-1]);

                fprintf(fout,"to ");

                for(k=0;k<Param.Dim-1;k++)
                    fprintf(fout,"%014.8f,\t",graph.r[j][k]);

                fprintf(fout,"%014.8f\n",graph.r[j][Param.Dim-1]);
            }
            fprintf(fout,"replot\n");
            fclose(fout);
}

void compute_degree(void)
{
    int il;

    for(il=0;il<graph.N_par;il++)
        graph.grado[il]=0;

    for(il=0;il<graph.N_links;il++)
        {
            graph.grado[graph.origen[il]]+=graph.peso[il]/2.0; //Para grado out
            graph.grado[graph.destino[il]]+=graph.peso[il]/2.0; //Para grado in
        }
}
