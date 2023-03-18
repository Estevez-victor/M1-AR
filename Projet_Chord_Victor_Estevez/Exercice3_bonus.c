#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define TAGINIT 1
#define TAGTAB 2
#define TAGINSERT 3
#define TAGCALC 4
#define TAGFIN 5


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
    for(int i=1; i<size-1; i++) {      // Tri du tableau
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

    int matFinger[size-1][M]; //matrice qui stockera les Finger Tables des pairs
    matFinger[0][0]=-1;

    int finger[M]; //table des finger
    int corr[M]; // correspondance Chord -> MPI
    for(int i=1; i<size; i++){
        MPI_Send(inUse, size-1, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); // envoye du tableau des pairs à i

        dpi=1;
        l=(i-1)%(size-1);
        for(int m=0; m<M; m++){      // calcule de la finger table de chaque processus
            if(m==0){
                tmp=(inUse[i-1]+1)%(max);
                finger[0]=inUse[i%(size-1)];
                corr[0]=(i%(size-1))+1;
            }else{
                next=-1;
                tmp=(tmp+dpi)%(max);
                dpi=dpi*2;
                
                while(next==-1){ // on tourne sur l'anneau jusqu'a trouver un intervale qui corespond a (pi+2^i)%2^M
                    if(l<size-2){ //dans un cas normal
                        if(tmp>inUse[l] && tmp<=inUse[l+1]){
                            next=inUse[l+1];  
                            corr[m]=l+2;
                        }
                    }else{ //dans le cas où on passe par le 0 de l'anneau
                        if(tmp>inUse[size-2] || tmp<=inUse[0]){
                            next=inUse[0];
                            corr[m]=1;
                        }
                    }
                    if(next==-1){
                        l=(l+1)%(size-1);
                    }
                }
                finger[m]=next;
            }  
        }
        printf("Finger table de %d:________\n",inUse[i-1]);
        for(int f=0;f<M;f++){
            printf("| %d(%d) ",finger[f],corr[f]);
        }
        printf("|\n");
        MPI_Send(&finger, M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); //envoye de la finger table au pair i
        MPI_Send(&corr, M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); // et leurs correspondances
        memcpy(matFinger[i], corr, 4*M);

    }



    int inverse[size-1];
    int pos;
    for(int i=1;i<size;i++){
        pos=0;
        for(int m=1;m<size;m++){
            for(int f=0;f<M;f++){
                if(matFinger[m][f]==i){
                    inverse[pos]=m;
                    pos++;
                    break;
                }
            }
        }
    MPI_Send(&pos, 1, MPI_INT, i, TAGINIT, MPI_COMM_WORLD);
    MPI_Send(&inverse, pos, MPI_INT, i, TAGINIT, MPI_COMM_WORLD);
    printf("inversep de %d [%d]:_____\n",i,pos);
        for(int i=0;i<pos;i++){
            printf("| %d ",inverse[i]);
        }
        printf("|\n");

    }
    MPI_Send(&(inUse[size-1]), 1, MPI_INT, size, TAGINIT, MPI_COMM_WORLD);// on envoye au nouveau son id chord,
    int up = 1 +( rand () % (nb_proc-2));// et le pair qu'il connaitera
    MPI_Send(&up, 1, MPI_INT, size, TAGINIT, MPI_COMM_WORLD);

}


void insertion(int nb_proc, int rang) {
    MPI_Status status;
    int p;
    int inversep[p];
    int inUse[nb_proc]; 
    int chord;
    int finger[M];        
    int corr[M];
    int up;
    int inUseCorr[nb_proc];

        int max=1;
        for(int i=0;i<M;i++){   // max = (2^M)
            max=max*2;
        }

    int run=1;
    if(rang<nb_proc-1){
        MPI_Recv(&inUse, nb_proc-2, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
        chord=inUse[rang-1];
        
        MPI_Recv(finger, M, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status); //reception de la Finger Table
        MPI_Recv(corr, M, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status); // et des correspondances

        MPI_Recv(&p, 1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
        MPI_Recv(&inversep, p, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);


    }else{
        MPI_Recv(&chord, 1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
        MPI_Recv(&up, 1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
        printf("\n[%d] je veux m'inserer.\n\n",chord);
        MPI_Send(&up, 1, MPI_INT, up, TAGTAB, MPI_COMM_WORLD);
        MPI_Recv(inUse, nb_proc-1, MPI_INT, up, TAGTAB, MPI_COMM_WORLD, &status);

   
        inUse[nb_proc-2]=chord; 
        for(int i=0;i<nb_proc;i++){
            inUseCorr[i]=i+1;
        }

        
        int j,k,c,d; //variables pour le tri
        for(int i=1; i<nb_proc-1; i++) {      // Tri du tableau
            if(inUse[i] < inUse[i-1]) {
                j = 0;
                while(inUse[j] < inUse[i]){
                    j++;
                }
                c = inUse[i];
                d = inUseCorr[i];
                for(k = i-1; k >= j; k--){
                    inUse[k+1] = inUse[k];
                    inUseCorr[k+1] = inUseCorr[k];
                }
                inUse[j] = c;
                inUseCorr[j] = d;
            }
        }
        
        int moi;
        for(int i=0;i<nb_proc-1;i++){
            if(inUse[i]==chord){
                moi=i;
            }
        }
        int next;
        int tmp;
        int dpi=1;
        int l=moi;
        for(int m=0; m<M; m++){      // calcule de la finger table de chaque processus
            if(m==0){
                tmp=(inUse[moi]+1)%(max);
                finger[0]=inUse[(moi+1)%(nb_proc-1)];
                corr[0]=((moi+1)%(nb_proc-1));
            }else{
                next=-1;
                tmp=(tmp+dpi)%(max);
                dpi=dpi*2;
                
                while(next==-1){ // on tourne sur l'anneau jusqu'a trouver un intervale qui corespond a (pi+2^i)%2^M
                    if(l<nb_proc-2){ //dans un cas normal
                        if(tmp>inUse[l] && tmp<=inUse[l+1]){
                            next=inUse[l+1];  
                            corr[m]=inUseCorr[l+1];
                        }
                    }else{ //dans le cas où on passe par le 0 de l'anneau
                        if(tmp>inUse[nb_proc-2] || tmp<=inUse[0]){
                            next=inUse[0];
                            corr[m]=inUseCorr[0];
                        }
                    }
                    if(next==-1){
                        l=(l+1)%(nb_proc-1);
                    }
                }
                finger[m]=next;
            }  
        }
        printf("Finger table de %d:________\n",inUse[moi]);
        for(int f=0;f<M;f++){
            printf("` %d(%d) ",finger[f],corr[f]);
        }
        printf("`\n");
        MPI_Send(&chord, 1, MPI_INT, corr[0], TAGINSERT, MPI_COMM_WORLD);
        run=0;
    }
    int mes;
    while(run){
        MPI_Recv(&mes, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(status.MPI_TAG==TAGTAB){
            MPI_Send(inUse, nb_proc-1, MPI_INT, status.MPI_SOURCE, TAGTAB, MPI_COMM_WORLD);
            continue;
        }
        if(status.MPI_TAG==TAGINSERT){
            int sent;
            for(int i=1;i<nb_proc-1;i++){
                for(int j=0;j<p;j++){
                    if(i==inversep[j]){
                        MPI_Send(&mes, 1, MPI_INT, i, TAGCALC, MPI_COMM_WORLD);
                        break;
                    }
                    if(j==p-1){
                        MPI_Send(&mes, 1, MPI_INT, i, TAGFIN, MPI_COMM_WORLD);
                    }
                }
            }
            run=0;
        }

        if(status.MPI_TAG==TAGCALC){
        
        inUse[nb_proc-2]=mes;
        
        int j,k,c,d; //variables pour le tri
        for(int i=1; i<nb_proc-1; i++) {      // Tri du tableau
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
        int moi;
        for(int i=0;i<nb_proc-1;i++){
            if(inUse[i]==chord){
                moi=i;
            }
        }
        int next;
        int tmp;
        int dpi=1;
        int l=moi;
        for(int m=0; m<M; m++){      // calcule de la finger table de chaque processus
            if(m==0){
                tmp=(inUse[moi]+1)%(max);
                finger[0]=inUse[(moi+1)%(nb_proc-1)];
            }else{
                next=-1;
                tmp=(tmp+dpi)%(max);
                dpi=dpi*2;
                
                while(next==-1){ // on tourne sur l'anneau jusqu'a trouver un intervale qui corespond a (pi+2^i)%2^M
                    if(l<nb_proc-2){ //dans un cas normal
                        if(tmp>inUse[l] && tmp<=inUse[l+1]){
                            next=inUse[l+1];  
                        }
                    }else{ //dans le cas où on passe par le 0 de l'anneau
                        if(tmp>inUse[nb_proc-2] || tmp<=inUse[0]){
                            next=inUse[0];
                        }
                    }
                    if(next==-1){
                        l=(l+1)%(nb_proc-1);
                    }
                }
                finger[m]=next;
            }  
        }
        printf("NOUVELLE Finger table de %d:________\n",chord);
        for(int f=0;f<M;f++){
            printf("` %d ",finger[f]);
        }printf("'\n");
        run=0;

        }
        if(status.MPI_TAG==TAGFIN){
            run=0;
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
        insertion(nb_proc, rang);
    }
    
    MPI_Finalize();
    return 0;
}
