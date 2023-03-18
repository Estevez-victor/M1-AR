#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TAGINIT 1
#define TAGCHERCHE 2
#define TAGTROUVE 3
#define TAGFIN 4

// IMPORTANT (K dans le td)
#define M 7


void simulateur(int nb_proc) { //Initialisation du systeme
    srand(time(NULL));
    int inUse[nb_proc-1];
    
    int max=1;
    for(int i=0;i<M;i++){   // max = (2^M)
        max=max*2;
    }
    printf("2^M: %d\n",max);

    int size=0;
    int randed=-1;
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

    int j,k,c; //variables pour le tri
    for(int i=1; i<size; i++) {      // Tri du tableau
        if(inUse[i] < inUse[i-1]) {
            j = 0;
            while(inUse[j] < inUse[i]){
                j++;
            }
            c = inUse[i];
            for(k = i-1; k >= j; k--){
                inUse[k+1] = inUse[k];
            }
            inUse[j] = c;
        }
    }
    printf("ID chord des pairs en jeu:\n");
    for(int i=0; i<size;i++){
        printf("| %d ",inUse[i]);
    }
    printf("|\n\n");
    int dpi; //utilisé pour le calcul de 2^i
    int tmp; //stock la valeur de (pi+2^i)%2^M
    int l; //position de recherche sur l'anneau
    int next; //utilisé pour arreter la boucle de recherche

    int finger[M]; //table des finger
    int corr[M]; // correspondance Chord -> MPI

    for(int i=1; i<=size; i++){
        MPI_Send(inUse, size, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); // envoye du tableau des pairs à i

        dpi=1;
        l=i-1;
        for(int m=0; m<M; m++){      // calcule de la finger table de chaque processus
            if(m==0){ //le premier element de la Finger Table est le successeur de i
                tmp=(inUse[i-1]+1)%(max);
                finger[0]=inUse[i%size];
                corr[0]=(i%size)+1;
            }else{
                next=-1;
                tmp=(tmp+dpi)%(max);
                dpi=dpi*2;
                
                while(next==-1){ // on tourne sur l'anneau jusqu'a trouver un intervale qui corespond a (pi+2^i)%2^M
                    if(l<size-1){ //dans un cas normal
                        if(tmp>inUse[l] && tmp<=inUse[l+1]){
                            next=inUse[l+1];
                            corr[m]=l+2;
                        }
                    }else{ //dans le cas où on passe par le 0 de l'anneau
                        if(tmp>inUse[size-1] || tmp<=inUse[0]){
                            next=inUse[0];
                            corr[m]=1;
                        }
                    }
                    if(next==-1){
                        l=(l+1)%size;
                    }
                }
                finger[m]=next;
            }  
        }
        printf("Finger table de %d:________\n",inUse[i-1]);
        for(int f=0;f<M;f++){
            printf("| %d ",finger[f]);
        }
        printf("|\n");
        MPI_Send(&finger, M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); //envoye de la finger table à i
        MPI_Send(&corr, M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); // et leurs correspondances
        
    }

    // On tire au hasard un chercheur et une clé et on lui envoye
    int atrouver = ( rand() % (max));
    int chercheur = 1 +( rand () % (nb_proc-1));
    int aps[2] = {atrouver, chercheur};   
    printf("\nA trouver:%d , Chercheur:%d  (note: -1 est le simulateur)\n\n",atrouver, inUse[chercheur-1]);
    MPI_Send(aps, 2, MPI_INT, chercheur, TAGCHERCHE, MPI_COMM_WORLD);
    
    MPI_Status status;
    MPI_Recv(aps, 2, MPI_INT, chercheur, TAGTROUVE, MPI_COMM_WORLD, &status); // On attend que le chercheur ait trouvé la clé
    printf("[Simulateur]Data bien recu de la part de %d, terminaison.\n", inUse[chercheur-1]);
    for(int i=1; i<nb_proc; i++){ //on envoye un message de terminaison a tout le monde
        MPI_Send(aps, 2, MPI_INT, i, TAGFIN, MPI_COMM_WORLD);
    }
}


void recherche(int nb_proc, int rang) {
    MPI_Status status;
    int inUse[nb_proc]; //decalage du tableau des id chord en jeu pour pouvoir printf facilement l'id chord des messages
    inUse[0]=-1;
    MPI_Recv(&inUse[1], nb_proc-1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
    int chord= inUse[rang];
    
    int finger[M];
    MPI_Recv(finger, M, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
    int corr[M]; // tableau de correspondance finger -> rang MPI
    MPI_Recv(corr, M, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);

    int elu=0; 
    int run=1;
    int aps[2]; // aps[valeur a chercher][processus MPI qui la cherche]
    int sent;
    while(run){
        MPI_Recv(aps, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(aps[1]==rang){// cas du processus chercheur
            elu=1;
        }
        if(status.MPI_TAG==TAGCHERCHE){
            printf("c[%d], <- demande de recherche recu de %d \n", chord, inUse[status.MPI_SOURCE]);
            sent=0;
            if(aps[0]==chord){ // Si notre id chord correspond a la data on en est responsable
                if(elu){
                    MPI_Send(aps, 2, MPI_INT, 0, TAGTROUVE, MPI_COMM_WORLD);
                    printf("c[%d], -> je suis responsable de la data, je l'envoye au simulateur.\n",chord);
                }else{
                    MPI_Send(aps, 2, MPI_INT, aps[1], TAGTROUVE, MPI_COMM_WORLD);
                    printf("c[%d], -> je suis responsable de la data, je la retourne à %d.\n", chord, inUse[aps[1]]);
                }
                continue;

            }else{ // On regarde si notre successeur est responsable de la data
                if(chord<finger[0]){// cas normal
                    if(aps[0]>chord && aps[0]<finger[0]){
                        MPI_Send(aps,2, MPI_INT, corr[0], TAGTROUVE, MPI_COMM_WORLD);
                        printf("c[%d], -> %d est responsable de la data, je le previens.\n", chord, finger[0]);
                        sent=1;
                    }
                }else{// cas où l'on passe par 0
                    if(aps[0]>chord || aps[0]<finger[0]){
                        MPI_Send(aps,2,MPI_INT, corr[0], TAGTROUVE, MPI_COMM_WORLD);
                        printf("c[%d], -> %d est responsable de la data, je le previens.\n",chord ,finger[0]);
                        sent=1;
                    }
                }
                if(!sent){ //si aucun des deux cas precedents, il faut chercher dans la finger table un intervale qui correspond
                    for(int i=0; i<M-1; i++){
                        if(finger[i]<=finger[i+1]){//cas normal
                            if(aps[0]>=finger[i] && aps[0]<finger[i+1]){
                                MPI_Send(aps, 2, MPI_INT, corr[i], TAGCHERCHE, MPI_COMM_WORLD);
                                printf("c[%d], en cours de recherche -> %d\n",chord, finger[i]);
                                sent=1;
                                break;
                            }
                        }else{//cas où l'on passe par 0
                            if(aps[0]>= finger[i] || aps[0]<finger[i+1]){
                                MPI_Send(aps, 2, MPI_INT, corr[i], TAGCHERCHE, MPI_COMM_WORLD);
                                printf("c[%d], en cours de recherche -> %d\n",chord, finger[i]);
                                sent=1;
                                break;
                            }
                        }
                    }
                }

                if(!sent){ // Et si aucun intervale dans la finger table ne correpond à notre data, on envoye au pair le plus lointain
                    MPI_Send(aps,2, MPI_INT, corr[M-1], TAGCHERCHE, MPI_COMM_WORLD);
                                printf("c[%d], en cours de recherche -> %d\n",chord, finger[M-1]);
                }
            }

            if(elu){ // le processus chercheur de met en attente la data et previendra le simulateur
                MPI_Recv(aps, 2, MPI_INT, MPI_ANY_SOURCE, TAGTROUVE, MPI_COMM_WORLD, &status);
                printf("c[%d], <- J'ai la data. \n",chord);
                MPI_Send(aps, 2, MPI_INT, 0, TAGTROUVE, MPI_COMM_WORLD);
                printf("c[%d], -> Je l'envoye au simulateur.\n", chord);
            }
            continue;
        }

        if(status.MPI_TAG==TAGTROUVE){ // on a été notifié par le precedent qu'on est responssable de la data
            MPI_Send(aps, 2, MPI_INT, aps[1], TAGTROUVE, MPI_COMM_WORLD);
            printf("c[%d], -> je suis responsable de la data, je la retourne à %d.\n", chord, inUse[aps[1]]);
            continue;
        }

        if(status.MPI_TAG==TAGFIN){ // message d'arret du simulateur
            run=0;
            continue;
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
        recherche(nb_proc, rang);
    }
    
    MPI_Finalize();
    return 0;
}
