
//
//  main.c
//  cammini
//
//  Created by Monica  on 04/05/23.
// CALCOLA, A VALORE FISSATO DI ETA, I CAMMINI A N VIA VIA DECRESCENTE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran2(void);
void ranstart(void);
void ranfinish(void);

void geometry(int, int*);
void initialize_lattice(int, int, float*);
void update_metropolis(int, int*, float*);
float measure(int nlatt,int *val, float* field);


long idum = 0, idum2 = 0, iy = 0, iv[NTAB];


    
// SUBROUTINES
    
    
// Definisco la cordinata + e -
    
void geometry(int nlatt, int* val){
  
    for(int i=0; i<nlatt-1; ++i){
        val[i]=i+1;
    }
    val[nlatt-1]=0;

    for(int i=nlatt+1; i<2*nlatt; ++i){
        val[i]=(i-1)-nlatt;
    }
    val[nlatt]=nlatt-1;
       
}
    
// Assegno la configurazione di partenza alla catena di Markow
    
void initialize_lattice(int iflag, int nlatt, float *field) {
  
    double x;
    
    FILE *lattice_file;
   
    // matrice con le variabili di spin
    
    lattice_file = fopen("/Users/monicacesario/Desktop/modulo3/path_int/lattice3.txt", "r");
    if (lattice_file == NULL) {
        printf("Errore: il file lattice non può essere aperto.\n");
        
    }
   
    //partenza a freddo
    if(iflag <= 0){
        for(int i = 0; i<nlatt; i++){
            
                field[i]=0.0;
            
        }
    }
    
    //partenza a caldo
    else if(iflag == 1){
        
        for(int i=0; i<nlatt; i++){
            
                //nro random tra -1 e 1
            x = 1.0 - 2.*(double)(rand()/(RAND_MAX + 1.0));
                field[i] = x;
                
            }
        
    }  // da dove ero rimasto l'ultima volta
    else {
        for (int i=0; i<nlatt; i++){
            
                fscanf(lattice_file, "%f", &field[i]);
            
        }
       
        }
       
     
        fclose(lattice_file);

    }
    
    
    
void update_metropolis(int nlatt, int *val, float *field){
       
        int npp[nlatt], nmm[nlatt];
        int ip,im;
        float force,phi,phi_prova;
        double p_rat, eta = 1.;
        float d_metro = 1.;
        double x,y;
       double c1 = 1./eta;
       double c2 = (1./eta + eta/2.);
        
        for(int i = 0; i < nlatt; i++){
            npp[i]=val[i];
        }

        for(int i=nlatt; i < 2*nlatt; i++){
            nmm[i-nlatt] = val[i];
        }
       
        //loop su tutti i siti
        for (int i = 0; i<nlatt; i++){
            
                ip = npp[i];
                im = nmm[i];
               
                
                force = field[ip]+field[im];
               
               // valore attuale del campo
                x = (double)(rand()/(RAND_MAX + 1.0));
                phi = field[i];
                phi_prova = phi + 2.*d_metro*(0.5-x);
                
               
                p_rat = c1 * phi_prova * force - c2 * pow(phi_prova,2);
            p_rat = p_rat - c1 * phi * force + c2 * pow(phi,2);
               
                // METRO-TEST y = random (0,1)
            y=log((double)(rand()/(RAND_MAX + 1.0)));
                
                // x<p_rat verifica anche il caso p_rat > 1, se sì accetto
                if(y<= p_rat){
                    field[i] = phi_prova;
                    
                }
                
                
            }
        }
       
        

    

    
float measure(int nlatt,int *val, float* field){
    
    float obs1 = 0.0;
    int npp[nlatt];
   
   // float obs2 = 0.0;
 
    
    
    
   // for(int i = 0; i < nlatt; i++){
     //   npp[i]=val[i];
    //}

   
    for(int i = 0; i < nlatt; i++){
        
        obs1 = obs1 + pow(field[i],2);
        //j = npp[i];
        //obs2 = obs2 + pow((field[i] - field[j]),2);
      
        
        
        
    }
    
   return obs1/(float)nlatt;
   
    
    
   
    
    
}


float ran2(void)

{

    int j;

    long k;

    float temp;

    if (idum <= 0) {

        if (-(idum) < 1) idum=1;

        else idum = -(idum);

        idum2=(idum);

        for (j=NTAB+7;j>=0;j--) {

            k=(idum)/IQ1;

            idum=IA1*(idum-k*IQ1)-k*IR1;

            if (idum < 0) idum += IM1;

            if (j < NTAB) iv[j] = idum;

        }

        iy=iv[0];

    }

    k=(idum)/IQ1;

    idum=IA1*(idum-k*IQ1)-k*IR1;

    if (idum < 0) idum += IM1;

    k=idum2/IQ2;

    idum2=IA2*(idum2-k*IQ2)-k*IR2;

    if (idum2 < 0) idum2 += IM2;

    j=iy/NDIV;

    iy=iv[j]-idum2;

    iv[j] = idum;

    if (iy < 1) iy += IMM1;

    if ((temp=AM*iy) > RNMX) return RNMX;

    else return temp;

}

void ranstart(void) {
    FILE *fp;
    
    int i;
    
    fp = fopen("/Users/monicacesario/Desktop/modulo3/path_int/randomseed3.txt", "r");
    if (fp == NULL) {
        fprintf(stderr, "Errore nell'apertura del file randomseed\n");
        exit(1);
    }
    
    fscanf(fp, "%li", &idum);
    //printf("idum=%li\n",idum);
    fscanf(fp, "%li", &idum2);
    //printf("idum2=%li\n",idum2);
    for (i = 0; i < 32; i++) {
        fscanf(fp, "%li", &iv[i]);
        //printf("iv=%li\n",iv[i]);
    }
    fscanf(fp, "%li", &iy);
    //printf("iy=%li\n",iy);
    
    if (idum >= 0) {
        idum = -idum - 1;
    }
    
    fclose(fp);
}


void ranfinish(void) {
    FILE *fp;
    
    
    fp = fopen("/Users/monicacesario/Desktop/modulo3/path_int/randomseed3.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Errore nell'apertura del file randomseed\n");
        exit(1);
    }
    
    fprintf(fp, "%li\n", idum);
    fprintf(fp, "%li\n", idum2);
    for (int i = 0; i < 32; i++) {
        fprintf(fp, "%li\n", iv[i]);
    }
    fprintf(fp, "%li\n", iy);
    
    fclose(fp);
}



    int main(void) {
        
        
       // int nlatt=20; //inizializzare solo nel main
        float d_metro, eta, obs1 = 0.0;
        int iflag, measures, i_decorrel, i_term;
        FILE *input_file, *lattice_file,*output_file;;
        time_t t;
        /* Intializes random number generator */
          srand((unsigned) time(&t));
        // apertura del file da cui leggere i parametri della simulazione
        input_file = fopen("/Users/monicacesario/Desktop/modulo3/path_int/input3.txt", "r");
        
        if (input_file == NULL) {
            printf("Errore: il file input non può essere aperto.\n");
            return 1;
        }
        
        printf("PARAMETERS\n");
        // partenza/caldo/freddo/precedente
        fscanf(input_file, "%d", &iflag);
        printf("IFLAG: %d\n", iflag);
        
        //numero di misure
        fscanf(input_file, "%d",&measures);
        printf("MEASURES: %d\n", measures);
        
        //updating tra una misura e l'altra
        fscanf(input_file, "%d", &i_decorrel);
        printf("DECORREL: %d\n", i_decorrel);
        
        //passi di termalizzazione
        fscanf(input_file, "%d", &i_term);
        printf("ITERM: %d\n", i_term);
        
        //parametro del metropolis (sarebbe il delta che mi da l'intervallo in cui scelogo y_prova?)
        fscanf(input_file, "%f", &d_metro);
        printf("D_METRO: %f\n", d_metro);
        
        //parametro eta = omega * a
       fscanf(input_file, "%f", &eta);
        printf("ETA: %f\n", eta);
        
        printf("\n\n");
        
        
        output_file = fopen("/Users/monicacesario/Desktop/modulo3/path_int/output.txt", "a");
        if (output_file == NULL) {
            printf("Errore: il file di output non può essere aperto.\n");
            return 1;
        }
        
        lattice_file = fopen("/Users/monicacesario/Desktop/modulo3/path_int/lattice3.txt", "w");
        if (lattice_file == NULL) {
            printf("Errore: il file lattice non può essere aperto.\n");
            return 1;
        }
        // inizializzo generatore di numeri random
      //  ranstart();
        
       // for (float eta = 0.6; eta > 0.01; eta -= 0.02) {
            
           
           
            int  count = 0;
            float sum_obs1 = 0.0;
            int nlatt = (int)round(1.0 / eta);
            int val[2*nlatt];
            float *field = (float *)malloc(nlatt * sizeof(float*));
            d_metro = 2.*sqrt(eta);
              printf("\n nlatt: %d NLATT*ETA: %d\n",nlatt, (int)round(nlatt*eta));
        // inizializzo configurazione inziale
        initialize_lattice(iflag, nlatt, (float*)field);
        
        // inizializzo condizioni al bordo
        geometry(nlatt,(int*)val);
        
        
        
        // TEMRALIZZAZIONE
     for (int i = 0; i < i_term; i++) {
            
         update_metropolis(nlatt, (int*)val, (float*)field);
        }
        
       
        
        // SESSIONE ALL'EQUILIBRIO CON MISURE
        for(int i=0; i<measures; i++){
            
            
            update_metropolis(nlatt, (int*)val, (float*)field);
            
             if((i+1) % i_decorrel == 0){
            
           obs1 = measure(nlatt,(int*) val, (float*) field);
            sum_obs1 += obs1;
             count += 1;
           
                 
            }
           
        }
    fprintf(output_file,"%f %f\n",eta,sum_obs1/count);
        printf("eta: %f Avg Obs: %f\n",eta,sum_obs1/count);
        //    printf("count: %d", count);
        
        //salvo la configurazione per poter eventualmente ripartire
       
        for (int i=0; i<nlatt; i++){
                   fprintf(lattice_file, "%f\n", field[i]);
               }
               
               free(field);
      //     }
           
           // Prendo ultimo stato del generatore random
         //  ranfinish();
        
        fclose(input_file);
        fclose(lattice_file);
        fclose(output_file);
    
        printf("\nDONE!\n");
        return 0;
    
}
