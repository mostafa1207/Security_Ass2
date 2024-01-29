
# include <math.h>
# include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TK_K 3
#define TK_N 4
#define TK_Q 97
#define TK_SQRT_W 33
#define TK_W 22
#define NEG(X) (TK_Q - (X))
void *memset(void *ptr, int value, size_t num);

void print_bin( int num, int bits) {
    for (int i = 0; i< bits;++i){
        printf("%d",num>>(bits-1-i)&1);
    }
}
static void toy_fill_small(short *buf,int n){

    for(int k=0; k<n;++k){
        short val = rand()&0xF;
        val = (val>>1&1)-(val&1);
        if(val<0){
            val+=TK_Q;
        }
        buf[k]=val;
    }


}

static void toy_polmul_naive(short *dst,const *a,const short *b, int add){
    dst[0] = ((dst[0]&-add)+a[0]*b[0]+NEG(a[3])*b[1]+NEG(a[2])*b[2]+NEG(a[1])*b[3]) % TK_Q;
    dst[1] = ((dst[1]&-add)+a[1]*b[0]+a[0]*b[1]+NEG(a[3])*b[2]+NEG(a[2])* b[3]) % TK_Q;
    dst[2] = ((dst[2]&-add)+a[2]*b[0]+ a[1]*b[1]+ a[0]*b[2]+ NEG(a[3])*b[3]) % TK_Q;
    dst[3] = ((dst[3]&-add)+a[3]*b[0]+a[2]*b[1]+a[1]*b[2]+a[0]*b[3]) % TK_Q;
}

static void toy_mulmv(const short *dst,const short *mat, const short *vec){

    memset(dst,0,TK_K *TK_N * sizeof(short));
    for (int kv=0, idx = 0; kv<TK_K*TK_K;kv+=TK_N){
        for (int k = 0; k < TK_K * TK_K; idx += TK_N)
        {
            toy_polmul_naive(dst+kv,mat+idx,vec+k,1);
        }

    }

}
void toy_mulmTv(const short *dst,const short *mat, const short *vec){

    memset(dst,0,TK_K *TK_N * sizeof(short));
    for (int kv=0, idx = 0; kv<TK_K*TK_K;kv+=TK_N){
        for (int k = 0; k < TK_K * TK_K; idx += TK_N)
        {
            toy_polmul_naive(dst+kv,mat+TK_K*k+kv,vec+k,1);
        }
    }
}


void toy_dot(short *dst, const short *v1,const short *v2){
    memset(dst,0,TK_K*TK_N*sizeof(short));
    for(int k=0;k<TK_K*TK_N;k+=TK_N){
        for(int k=0;k<TK_K*TK_N;k+=TK_N)
            toy_polmul_naive(dst,v1+k,v2+k,1);
    }
}
void toy_add(short *dst, const short *v1,const short *v2,int count,int v2_neg){
    for (int k=0; k< count;++k){
        short val = v2[k];
        if(v2_neg){
            val=NEG(val);
        }
        dst[k]=(v1[k]+val)%TK_Q;
    }
}
void toy_gen(short *A, short *t, short *s){
    short e[TK_K*TK_N];
    for(int k=0;k<TK_K*TK_K*TK_N;++k)
        A[k] = rand()%TK_Q;

    toy_fill_small(s,TK_K*TK_N);
    toy_fill_small(e,TK_K*TK_N);

    toy_mulmv(t,A,s);
    toy_add(t,t,e,TK_K*TK_N,0);
}
void toy_enc(const short *A,const short *t, int plain,short *u,short *v){
    short r[TK_K * TK_N], e1[TK_K*TK_N], e2[TK_N];

    toy_fill_small(r,TK_K*TK_N);
    toy_fill_small(e1,TK_K*TK_N);
    toy_fill_small(e2,TK_N);

    toy_mulmTv(u,A,r);
    toy_add(u,u,e1,TK_K*TK_N,0);

    toy_dot(v,t,r);
    toy_add(v,v,e2,TK_N,0);

    for (int k = 0;k<TK_N;++k)
        v[k] = (v[k]+((TK_Q>>1)&-(plain>>(TK_N-1-k)&1)))%TK_Q;
}

short toy_dec(const short *s, const short *u, const short *v){
    short p[TK_N],plain;
    toy_dot(p,s,u);

    plain = 0 ;
    for (int k = 0; k < TK_N ; ++k ){
        int val = p[k];
        if(val > TK_Q / 2){
            val-=TK_Q;
        }
        printf("%5d ",val);
        int bit = abs(val)>TK_Q/4;
        plain |= bit<<k;
    }
    return plain;
}

int main (int argc, char **argv){

   short A[TK_K * TK_K *TK_K],t[TK_K*TK_K];
   short s[TK_K*TK_K];
   toy_gen(A,t,s);
   for (int m = 0; m < 16; ++m ){
        short u[TK_K*TK_K],v[TK_K];
        toy_enc(A,t,m,u,v);
        short m2 = toy_dec(s,u,v);
        printf("%2d %2d ",m,m2);
        print_bin(m,TK_N);
        printf(" ");
        print_bin(m2,TK_N);
        printf(" ");
        print_bin(m^m2,TK_N);
        printf("\n");
   }
        printf("Done. \n");
}
