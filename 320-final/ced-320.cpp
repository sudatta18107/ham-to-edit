// Generates a collection of strings satisfying property 3a

#include <sys/time.h>  
#include <stdio.h>  
#include <string.h>
#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <algorithm>
#include <iostream>
 

#define MAXLEN 100000
#define K  (MAXLEN)

typedef int64_t int64;
typedef int32_t int32;
typedef int16_t int16;
typedef int8_t int8;
#define max(x,y)  (((x)>(y))?(x):(y))
#define min(x,y)  (((x)<(y))?(x):(y))


#define LEN  320
#define RATE 8
#define COLSIZE 1500
#define RANKFRAC 0.0001

int64  collection[COLSIZE];             // compressed collection strings
char  collectionstr[COLSIZE][LEN+1];    // C-like strings from the collection (uncompressed)
char  collectionstrr[COLSIZE][LEN+1];   // reverse of all the strings from the collection

int64  collfail[COLSIZE];   // failure statistics

int cnt=0;     // current size of the collection

int Pstat[LEN+1][LEN+1];
int Mstat[(LEN+3)*(LEN+2)/2][LEN+1];
int Sstat[LEN+1][LEN+1];

int Pthr[LEN+1];
int Mthr[(LEN+3)*(LEN+2)/2];
int Sthr[LEN+1];



#include "Randombits.h"
uint32_t rndxy, rndprev, rndcnt;

int32 random32()
{
	rndxy++;
        if((rndcnt+=3)<3)rndcnt++;
	rndprev ^= (1|*((int32*)(rndtbl+(rndcnt&0xFFFF)))) * ((*((int32*)(rndtbl+(rndxy & 0xFFFF)))) ^ (*((int32*)(rndtbl+((rndxy >> 16)&0xFFFF)))) ^ (*((int32*)(rndtbl+((rndcnt >> 16)&0xFFFF)))) ^ (*((int32*)(rndtbl+((rndprev>>16) & 0xFFFF)))));
	return rndprev;
}

int64 random64()
{
    return ((int64)random32()) ^ (((int64)random32())<<32) ^ (((int64)random32())<<16); 
}


void srandom32(int32 s=0)
{
	rndcnt=s*13564;
	rndprev=10;
	rndxy=s;
}



int ED(const char *x, const int x_len, const char *y, const int y_len)
{
    if(x_len > y_len) return ED(y, y_len, x, x_len);
    if(x_len==0) return y_len;

    int fc_buf[K+2], fp_buf[K+2];          
    int *fc = fc_buf + 1;
    int *fp = fp_buf + 1;                  
    int i, j;

    
    for(j = 0; j <= x_len; j++) {
        fp[j] = j;
        fc[j] = j;
    }

    for(i = 1; i <= y_len; i++) {
        std::swap(fc, fp);
        fc[0] = i;

        for(j = 1; j <= x_len; j++) {
            if(x[j-1] == y[i-1]) {
                //Match
                fc[j] = fp[j - 1];
            } else {
                // edit operations
                fc[j] = 1 + min(min(fp[j], fc[j - 1]), fp[j - 1]);
            }
        }
    }

    return fc[x_len];
}

int EDstarmin(const char *y, const int y_len, const char *x, const int x_len)
// returns the edit distance of the best prefix match of x to y.
{
    int fc_buf[K+2], fp_buf[K+2];          
    int *fc = fc_buf + 1;
    int *fp = fp_buf + 1;                  
    int i, j;

    
    for(j = 0; j <= x_len; j++) {
        fp[j] = j;
        fc[j] = j;
    }

    for(i = 1; i <= y_len; i++) {
        std::swap(fc, fp);
        fc[0] = i;

        for(j = 1; j <= x_len; j++) {
            if((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*')) {
                //Match
                fc[j] = fp[j - 1];
            } else {
                // edit operations
                fc[j] = 1 + min(min(fp[j], fc[j - 1]), fp[j - 1]);
            }
        }
    }

    int res=x_len+y_len;
    for(j=0;j<=x_len;j++)res = min(res, fc[j]); 

    return res;
}


int EDstarpref(const char *y, const int y_len, const char *x, const int x_len, int *res)
// returns the edit distance of all prefixes of y to x.
{

    int fc_buf[K+2], fp_buf[K+2];          
    int *fc = fc_buf + 1;
    int *fp = fp_buf + 1;                  
    int i, j;

    
    for(j = 0; j <= x_len; j++) {
        fp[j] = j;
        fc[j] = j;
    }

    res[0]=x_len;

    for(i = 1; i <= y_len; i++) {
        std::swap(fc, fp);
        fc[0] = i;

        for(j = 1; j <= x_len; j++) {
            if((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*')) {
                //Match
                fc[j] = fp[j - 1];
            } else {
                // edit operations
                fc[j] = 1 + min(min(fp[j], fc[j - 1]), fp[j - 1]);
            }
        }

        res[i]=fc[x_len];
    }

    return fc[x_len];
}


int EDstar(const char *x, const int x_len, const char *y, const int y_len)
{
 //   if(x_len > y_len) return ED(y, y_len, x, x_len);
    if(x_len==0) return y_len;

    int fc_buf[K+2], fp_buf[K+2];          
    int *fc = fc_buf + 1;
    int *fp = fp_buf + 1;                  
    int i, j;

    
    for(j = 0; j <= x_len; j++) {
        fp[j] = j;
        fc[j] = j;
    }

    for(i = 1; i <= y_len; i++) {
        std::swap(fc, fp);
        fc[0] = i;

        for(j = 1; j <= x_len; j++) {
            if((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*')) {
                //Match
                fc[j] = fp[j - 1];
            } else {
                // edit operations
                fc[j] = 1 + min(min(fp[j], fc[j - 1]), fp[j - 1]);
            }
        }
    }

    return fc[x_len];
}



void randomstring(char *s, int len)
{
    for(;len>0;s++,len--)if((random32()&0x800)==0)*s='0';else *s='1';

    *s='\0';
}

void biasstr(char *s, int t,char a, char b)
{
   int i,j;
  
   for(i=j=0;i<t;i++)if(s[i]==a)j++;

   if(j<=1)return;

   j=(random32()&&0xfffffff)%j;

   for(i=0;i<t;i++)if(s[i]==a){
      if(j--==0){s[i]=b;return;}
   }
}


void inttostr(int64 x,char *s, int len)
{
    for(int i=0;i<len;s++,i++,x>>=1)if((x&0x01)==0)*s='0';else *s='1';
    *s = 0;
}

void strstar(char *s, int len, int t)
{
    for(int i=0;i<len;s++,i++)if((i%t) == 0)*s='*';
}

void reversestr(char *s, int len, char *r)
{
    for(int i=0;i<len;i++)r[len-1-i]=s[i];
    r[len]=0;
}

void randombiasstarstring(char *s)
{
    randomstring(s,LEN);
    strstar(s,LEN,RATE);
//    biasstr(s+1,RATE-1,'0','1');
//    biasstr(s+1,RATE-1,'0','1');
//    biasstr(s+RATE+1,RATE-1,'0','1');

//    biasstr(s+LEN-RATE,RATE-1,'1','0');
//    biasstr(s+LEN-RATE,RATE-1,'1','0');
//    biasstr(s+LEN-2*RATE,RATE-1,'1','0');
}


///////////////////////////////////////////////////////////////////////////////////////////////


void resetstat(int *stat)
{
   for(int i=0;i<=LEN;i++)stat[i]=0;
}

int rankstat(int *stat, int rank)
{
   int sum=0;

   for(int i=0;i<=LEN;i++)if((sum+=stat[i])>=rank)return i;
   return LEN;
}

void printstat(int *stat)
{
   for(int i=0;i<=LEN;i++)printf("%d ",stat[i]);
     printf("\n");
}

void stringstat(int64 cnt)
{
    int i,j,k,l,r;

    char x[LEN+1],y[LEN+1];
    char xr[LEN+1],yr[LEN+1];

    int res[LEN+1];

    for(i=0;i<=LEN;i++)resetstat(Pstat[i]);
    for(i=0;i<(LEN+3)*(LEN+2)/2;i++)resetstat(Mstat[i]);
    for(i=0;i<=LEN;i++)resetstat(Sstat[i]);

    j=0;

    while(j<cnt){

        randombiasstarstring(x);
        reversestr(x,LEN,xr);
        randombiasstarstring(y);
        reversestr(y,LEN,yr);

        for(l=0;l<LEN;l++){
           EDstarpref(x+l,LEN-l,y,LEN,res);
           for(r=l+1;r <= LEN; r++)Mstat[r*(r+1)/2+l][res[r]]++;
        }

       for(l=0;l<LEN;l++)Pstat[l][EDstarmin(xr,LEN-l,yr,LEN)]++;
    
       for(l=0;l<LEN;l++)Sstat[l][EDstarmin(x+l,LEN-l,y,LEN)]++;
    
      j++;
      if((j%100)==0){printf("+");fflush(stdout);}
   }

//   printstat(Pstat[33]);
//  printstat(Pstat[LEN-33]);

//   for(l=0;l<LEN;l++){ printf("* %d  %d  %d\n",rankstat(Pstat[l],RANKFRAC*cnt),rankstat(Sstat[LEN-l],RANKFRAC*cnt),rankstat(Pstat[l],RANKFRAC*cnt)+rankstat(Sstat[LEN-l],RANKFRAC*cnt)); /* printstat(Sstat[l]); */}
       
     for(l=0;l<LEN;l++){ Pthr[l]=rankstat(Pstat[l],RANKFRAC*cnt); Sthr[l]=rankstat(Sstat[l],RANKFRAC*cnt); }

     for(l=0;l<LEN;l++)for(r=l+1;r <= LEN; r++)Mthr[r*(r+1)/2+l] = rankstat(Mstat[r*(r+1)/2+l],RANKFRAC*cnt);

     int mined=LEN;

     for(l=2;l<LEN;l++)
        for(r=l+1;r<LEN;r++)mined=min(mined, Pthr[LEN-l] + Mthr[r*(r+1)/2+l] + Sthr[r]);
 
     for(l=0;l<LEN;l++)mined = min(mined, Pthr[l] + Sthr[LEN-l]);

     printf("\nmin = %d\n",mined);

     printf("int Pthr[]={");for(l=0;l<LEN;l++)printf("%d,",Pthr[l]);printf("-1};\n");
     printf("int Sthr[]={");for(l=0;l<LEN;l++)printf("%d,",Sthr[l]);printf("-1};\n");
     printf("int Mthr[]={");for(l=0;l<(LEN+3)*(LEN+2)/2;l++)printf("%d,",Mthr[l]);printf("-1};\n");
     

}


int test(char *x, char *y)
{
    int i,j,k,l,r;

    char xr[LEN+1],yr[LEN+1];

    int res[LEN+1];

    reversestr(x,LEN,xr);
    reversestr(y,LEN,yr);

    for(l=0;l<LEN;l++){
       EDstarpref(x+l,LEN-l,y,LEN,res);
       for(r=l+1;r <= LEN; r++)if(Mthr[r*(r+1)/2+l] > res[r])return -1;
    }

    for(l=0;l<LEN;l++)if(Pthr[l]> EDstarmin(xr,LEN-l,yr,LEN))return -1;
    
    for(l=0;l<LEN;l++)if(Sthr[l] > EDstarmin(x+l,LEN-l,y,LEN)) return -1;
    
    return 1;
}

int main(int argc, char* argv[])
{

   int i,j,k,l,r,offs;
   int testres;
   int rep=0;
   int remcnt;

   printf("LEN=%d RATE=%d RANKFRAC=%f\n", LEN, RATE, RANKFRAC);

   if(argc>1){
      int seed=12437*std::stoi(argv[1]);
      srandom32(seed);
      printf("Random seed=%d\n",seed);
   }
 
   stringstat(300000);
   
   for(remcnt=offs=cnt=0;cnt<COLSIZE;){
      randombiasstarstring(collectionstr[cnt]);
      
      testres=1;
      for(j=0;(j<cnt)&&(testres>0);j++)testres=min(test(collectionstr[(j+offs)%cnt],collectionstr[cnt]),test(collectionstr[cnt],collectionstr[(j+offs)%cnt]));

      if(testres==1){
         printf("\n%4d %s\n",cnt,collectionstr[cnt]);fflush(stdout);
         cnt++;rep=0;collfail[cnt]=0;remcnt++;

         int maxfail=0;
         int maxk=0;
         for(k=0;k<cnt;k++){ printf("%ld ",collfail[k]); if(maxfail<collfail[k]){maxk=k; maxfail=collfail[k];} } printf("\n");
         if(remcnt>=10){
           printf("removing %d\n",maxk);
           strcpy(collectionstr[maxk],collectionstr[cnt-1]);
           collfail[maxk]=0;
           cnt--;remcnt=0;
           for(k=0;k<cnt;k++)printf("col[%4d]=\"%s\";\n",k,collectionstr[k]);fflush(stdout);
         }
      }else{
        collfail[(j+offs)%cnt]++;
        offs=(1+j+offs)%cnt;
        if((rep++%100)==0){printf("+");fflush(stdout);}
      }
  }

   return 0;
}


