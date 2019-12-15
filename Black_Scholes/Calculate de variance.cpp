
int MeanT(int T[],int start,int finish){
    int m=0;
    for(int i=start,i<finish,i++){
        m+=T[i];
        }
    m=(m/(finish-start));
    return(m)
}
;


int ecartT(int T[],int start,int finish){
    int s=0;
    for(int i=start,i<finish,i++){
        s+=((T[i]-m)*2(T[i]-m));
    };
    s=(s/(finish-start));
    return(s)
}
;
