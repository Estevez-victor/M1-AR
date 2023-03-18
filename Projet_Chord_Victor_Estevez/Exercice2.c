#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TAGINIT 1
#define TAG_OUT 2
#define TAG_IN 3
#define TAGFIN 4
#define TAGTAB 5

// IMPORTANT (K dans le td)
#define M 8


void simulateur(int nb_proc) {
    srand(time(NULL));
    int inUse[nb_proc-1];
    int randed, size, max;
    max=1;
    for(int i=0;i<M;i++){   // max = (2^M)
        max=max*2;
    }
    printf("2^M: %d\n",max);
    
    size=0; 
    randed=-1;
    while(size<nb_proc-1){  // On tire aléatoirement une id chord unique pour chaque pair
        while(randed==-1){
            randed = ( rand() % (max));
            for(int i=0; i<size;i++){
                if(randed==inUse[i]){
                    randed=-1;
                    i=size;
                }
            }
            if(randed!=-1){
                inUse[size]=randed;
            }
        }
        size++;
        randed=-1;
    }
    
    printf("ID chord des pairs en jeu:\n");
    for(int i=0; i<size;i++){
        printf("| %d ",inUse[i]);
    }
    printf("|\n\n");
    
    randed=0;
    int initTab[size];
    for(int i=0; i<size; i++){ //On tire au hasard qui est initiateur
        initTab[i]=( rand() % 4 ); // un chance sur 4 d'etre initiateur
        if(!(initTab[i])){
            initTab[i]=1;
            randed=1;
        }else{
            initTab[i]=0;
        }
    }
    if(!randed){// on s'assure qu'il y ait au moins un initiateur
        initTab[0]=1;
    }
    int message[2];
    for(int i=1; i<=size; i++){
        message[0]=inUse[i-1];
        message[1]=initTab[i-1];
        MPI_Send(message, 2, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); // envoye de son id CHORD attitré et son status d'initiateur
    }

    
}


void calculFinger(int nb_proc, int rang) {
    MPI_Status status;
    int message[2];
    MPI_Recv(message, 2, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
    int chord=message[0];
    int initiateur= message[1];

    int suivant= (rang+1)%nb_proc;
    if(!suivant){
        suivant=1;
    }
    int precedent= (rang -1);
    if(!precedent){
        precedent=nb_proc-1;
    }
    
    int run=1;
    int distance=1; //la distance que le message devra parcourir
    int a;        //sauvegarde temporaire la distance du tour actuel
    int backIn=2; //on doit attendre les deux in pour pouvoir passer au tour suivant proprement
    int jeton[2]; // le jeton contient: {distance restante a parcouris, rang du processus qui a envoyé le jeton};
    int elu=0; // l'id MPI du processus élu

    while(run){   // algorthme d'election d'un leader: [Hirschberg & Sinclair]
        if(initiateur && (backIn==2)){
            jeton[0]=distance;
            jeton[1]=rang;
            MPI_Send(jeton, 2, MPI_INT, suivant, TAG_OUT, MPI_COMM_WORLD);
            MPI_Send(jeton, 2, MPI_INT, precedent, TAG_OUT, MPI_COMM_WORLD);
            a=distance;
            distance=distance*2;
            backIn=0;
        }
        MPI_Recv(jeton, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(status.MPI_TAG==TAG_OUT){
            if(initiateur){
                if(jeton[1]>rang){
                    initiateur=0;
                }else{
                    if(jeton[1]==rang){
                        for(int i=1; i<nb_proc; i++){
                            MPI_Send(jeton, 2, MPI_INT, i, TAGFIN, MPI_COMM_WORLD);
                        }
                    }
                }
            }
            if(!initiateur){
                jeton[0]=jeton[0]-1;
                if(jeton[0]>0){
                    if(status.MPI_SOURCE==precedent){
                        MPI_Send(jeton, 2, MPI_INT, suivant, TAG_OUT, MPI_COMM_WORLD);
                    }else{
                        MPI_Send(jeton, 2, MPI_INT, precedent, TAG_OUT, MPI_COMM_WORLD);
                    }
                }else{
                    if(status.MPI_SOURCE==precedent){
                        MPI_Send(jeton, 2, MPI_INT, precedent, TAG_IN, MPI_COMM_WORLD);
                    }else{
                        MPI_Send(jeton, 2, MPI_INT, suivant, TAG_IN, MPI_COMM_WORLD);
                    }
                }
            }
        }
        if(status.MPI_TAG==TAG_IN){
            if(initiateur){
                backIn=backIn+1;
            }else{
                if(status.MPI_SOURCE==precedent){
                        MPI_Send(jeton, 2, MPI_INT, suivant, TAG_IN, MPI_COMM_WORLD);
                }else{
                        MPI_Send(jeton, 2, MPI_INT, precedent, TAG_IN, MPI_COMM_WORLD);
                }
            }
        }
        if(status.MPI_TAG==TAGFIN){
            elu=status.MPI_SOURCE;
            run=0;
        }
    }//A la fin de l'algorithme un processus a été élu, et tout le mondes sait de qui il s'agit

    int finger[M];   
    if(elu!=rang){ //les processus lambda
        MPI_Send(&chord, 1, MPI_INT, elu, TAGTAB, MPI_COMM_WORLD); // chaque processus envoye son id chord a l'elu
        MPI_Recv(finger, M, MPI_INT, elu, TAGTAB, MPI_COMM_WORLD, &status); //puis attend sa table.
        printf("je suis %d, ma Finger Table:\n",chord);
        for(int i=0;i<M;i++){
            printf("| %d ",finger[i]);
        }
        printf("|\n\n");
    }else{ // le processus élu
        int chordTab[nb_proc-1]; 
        int corr[nb_proc-1];
        chordTab[rang-1]=chord;
        corr[nb_proc-2]=nb_proc-1;
        for(int i=1; i<nb_proc-1; i++){
            MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE,TAGTAB, MPI_COMM_WORLD, &status);
            chordTab[status.MPI_SOURCE-1]=a;
            corr[i-1]=i;
        }

        int j,k,c,d; //variables pour le tri
        for(int i=1; i<nb_proc-1; i++) {      // Tri du tableau
            if(chordTab[i] < chordTab[i-1]) {
                j = 0;
                while(chordTab[j] < chordTab[i]){
                    j++;
                }
                c = chordTab[i];
                d = corr[i];
                for(k = i-1; k >= j; k--){
                    chordTab[k+1] = chordTab[k];
                    corr[k+1] = corr[k];
                }
                chordTab[j] = c;
                corr[j] = d;
            }
        }
        printf("ID chord trié par le processus élu(%d):\n",rang);
        for(int i=0; i<nb_proc-1;i++){
            printf("| %d ",chordTab[i]);
        }
        printf("|\n\n");
        int max=1;
        for(int i=0;i<M;i++){   // max = (2^M)
            max=max*2;
        }
        int dpi; //utilisé pour le calcul de 2^i
        int tmp; //stock la valeur de (pi+2^i)%2^M
        int l; //position de recherche sur l'anneau
        int next; //utilisé pour arreter la boucle de recherche
        int size=nb_proc-1;

        for(int i=1; i<=size; i++){

            dpi=1;
            l=(i-1)%size;
            for(int m=0; m<M; m++){      // calcule de la finger table de chaque processus
                if(m==0){ //le premier element de la Finger Table est le successeur de i
                    tmp=(chordTab[i-1]+1)%(max);
                    finger[0]=chordTab[i%size];
                }else{
                    next=-1;
                    tmp=(tmp+dpi)%(max);
                    dpi=dpi*2;
                    
                    while(next==-1){
                        if(l<size-1){ //dans un cas normal
                            if(tmp>chordTab[l] && tmp<=chordTab[l+1]){
                                next=chordTab[l+1];
                            }
                        }else{ //dans le cas où on passe par le 0 de l'anneau
                            if(tmp>chordTab[size-1] || tmp<=chordTab[0]){
                                next=chordTab[0];
                            }
                        }
                        if(next==-1){
                            l=(l+1)%size;
                        }
                    }
                    finger[m]=next;
                }  
            }
            MPI_Send(&finger, M, MPI_INT, corr[i-1], TAGTAB, MPI_COMM_WORLD); //envoye de la finger table
        }


    }

}



/******************************************************************************/

int main (int argc, char* argv[]) {
    int rang, nb_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);

    
    if (rang == 0) {
        printf("M=%d\n",M);
        simulateur(nb_proc);
    } else {
        calculFinger(nb_proc, rang);
    }
    
    MPI_Finalize();
    return 0;
}
