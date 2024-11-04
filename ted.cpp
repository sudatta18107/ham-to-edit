// Generates data

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


#define LEN  54
#define RATE 18
#define COLSIZE 1200
#define RANKFRAC 0.001
#define ALPHA  0.14

int64  collection[COLSIZE];             // compressed collection strings
char  collectionstr[COLSIZE][LEN+1];    // C-like strings from the collection (uncompressed)
char  collectionstrr[COLSIZE][LEN+1];   // reverse of all the strings from the collection

int64  collfail[COLSIZE];   // failure statistics

int cnt=0;     // current size of the collection

#define XCOLSIZE 10
int8 LP[XCOLSIZE][XCOLSIZE][LEN+1];
int8 LM[XCOLSIZE][XCOLSIZE][(LEN+3)*(LEN+2)/2];
int8 LS[XCOLSIZE][XCOLSIZE][LEN+1];

int Pstat[LEN+1][LEN+1];
int Mstat[(LEN+3)*(LEN+2)/2][LEN+1];
int Sstat[LEN+1][LEN+1];

int Pthr[LEN+1];
int Mthr[(LEN+3)*(LEN+2)/2];
int Sthr[LEN+1];



#include "Randombits.h"
int32 rndxy, rndprev, rndcnt;

int32 random32()
{
	if(rndcnt--==0){rndcnt=105542531;rndprev++;}
	rndxy++;
	rndprev ^= (*((int32*)(rndtbl+(rndxy & 0xFFFF)))) ^ (*((int32*)(rndtbl+((rndxy >> 16)&0xFFFF)))) ^ (*((int32*)(rndtbl+(rndprev & 0xFFFF))));
	return rndprev;
}

int64 random64()
{
    return ((int64)random32()) ^ (((int64)random32())<<32) ^ (((int64)random32())<<16); 
}


void srandom32(int32 s=0)
{
	rndcnt=105542531;
	rndprev=10;
	rndxy=s;
}



int LCS(const char *x, const int x_len, const  char *y, const int y_len)
// computes LCS of x and y
{
	if(x_len > y_len)return LCS(y,y_len,x,x_len);

	int fc_buf[K+2],fp_buf[K+2];          // LCS buffers
	int *fc=fc_buf+1;
	int *fp=fp_buf+1;                     // current D and previous D
	int i,j;

	for(i=0;i<K+2;i++)fc_buf[i]=fp_buf[i]=0;

        for (i = 0; i < y_len; i++) {
           std::swap(fc,fp);
           
           for (j = 0; j < x_len; j++) {
               fc[j] = max(fp[j], fc[j - 1]);
               if (x[j] == y[i])fc[j] = max(fc[j],fp[j-1] + 1);
               
           }
       }

    return fc[x_len-1];
}

int slantedLCS(const char *x, const int x_len, const  char *y, const int y_len)
// computes LCS of x and y
{
	if(x_len > y_len)return LCS(y,y_len,x,x_len);

	int fc_buf[K+2],fp_buf[K+2];          // LCS buffers
	int *fc=fc_buf+1;
	int *fp=fp_buf+1;                     // current D and previous D
	int i,j;

	for(i=0;i<K+2;i++)fc_buf[i]=fp_buf[i]=0;

        for (i = 0; i < y_len; i++) {
           std::swap(fc,fp);
           
           for (j = 0; j < x_len; j++) {
               fc[j] = max(fp[j], fc[j - 1]);
               if ((x[j] == y[i])&& (i<j))fc[j] = max(fc[j],fp[j-1] + 1);
           }
       }

    return fc[x_len-1];
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


int EDstarprefnonstraight(const char *y, const int y_len, const char *x, const int x_len, int *res)
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
            if((i!=j)&&((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*'))) {
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


int EDstardiagprefnonstraight(const char *y, const int y_len, const char *x, const int x_len, int *res)
// returns the edit distance of all prefixes of y to prefixes of the same size of x.
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
            if((i!=j)&&((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*'))) {
                //Match
                fc[j] = fp[j - 1];
            } else {
                // edit operations
                fc[j] = 1 + min(min(fp[j], fc[j - 1]), fp[j - 1]);
            }
        }

        res[i]=fc[i];
    }

    return fc[x_len];
}


int EDstarnonstraight(const char *x, const int x_len, const char *y, const int y_len)
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
            if((i!=j)&&((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*'))) {
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

int EDstarminrangenonstraight(const char *y, const int y_len, const char *x, const int x_len, int x_low, int x_high, int x_diff)
// returns the edit distance of the best prefix match of y to x.
{
    int fc_buf[K+2], fp_buf[K+2];          
    int *fc = fc_buf + 1;
    int *fp = fp_buf + 1;                  
    int i, j;
    int res=x_len+y_len;

    
    for(j = 0; j <= x_len; j++) {
        fp[j] = j;
        fc[j] = j;
    }

    for(i = 1; i <= y_len; i++) {
        std::swap(fc, fp);
        fc[0] = i;

        for(j = 1; j <= x_len; j++) {
            if( ( (j<x_low)||(j>=x_high)||(j+x_diff!=i))  &&  ((x[j-1] == y[i-1])||(x[j-1]=='*')||(y[i-1]=='*'))) {
                //Match
                fc[j] = fp[j - 1];
            } else {
                // edit operations
                fc[j] = 1 + min(min(fp[j], fc[j - 1]), fp[j - 1]);
            }
        }

        res = min(res, fc[x_len]); 
    }


    return res;
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

const char *colstr[]={
/*   0*/  "*01000110110101100*00000000001000100*10010111011101001",
/*   1*/  "*01111100000010011*10010101110100101*01111110101000111",
/*   2*/  "*00111011010010000*11111101001100010*01101110101110110",
/*   3*/  "*00111100100111011*00000101110000111*10101111010100000",
/*   4*/  "*11101000010000001*11111101000011110*10000100001001000",
/*   5*/  "*11111110101101001*00001100011001001*00011001010000110",
/*   6*/  "*10000000001001000*10000001000111110*10011010100010100",
/*   7*/  "*10000001001001110*10111101111101000*11100000100101101",
/*   8*/  "*10111111101110101*00101010101110001*10001011000101110",
/*   9*/  "*01110100001100110*11010011110010100*00110010011110000",
/*  10*/  "*10101110111100011*11011101011111111*01000110100100001",
/*  11*/  "*11100011011110101*00001101000110011*01111011100001101",
/*  12*/  "*11111101111001101*00001101000011111*11011111111100111",
/*  13*/  "*11000100011000100*11110111111001110*00000000011001001",
/*  14*/  "*11001000001011110*01110000011101111*00011010000000001",
/*  15*/  "*00110101111110011*10100111111110110*11000001111100101",
/*  16*/  "*10101111100010011*01110110110110000*01110010000010011",
/*  17*/  "*10001000000010101*01000011111101110*10011010011011111",
/*  18*/  "*00010000100011101*00000000100100000*10111011011110100",
/*  19*/  "*01100111100011111*00100010100000000*01000110110011000",
/*  20*/  "*10110100001000111*10000111011110000*10001111001111011",
/*  21*/  "*10000101001011011*10110011111000110*10111011000010001",
/*  22*/  "*01011111110101111*11111000001111010*10111111000111110",
/*  23*/  "*00101101110000001*01010111011110000*11010101001110001",
/*  24*/  "*11111100101100010*10101010011101000*01000000011101101",
/*  25*/  "*11010000101100001*01111110000000001*11111110000001110",
/*  26*/  "*00111011111010000*00100111101101110*11000111111101100",
/*  27*/  "*11110010010000010*01100111100111111*10000001111001100",
/*  28*/  "*01010110011001111*11011110010111101*01000000000010111",
/*  29*/  "*00111011111111101*10111101100110010*01111000000000000",
// /*  30*/  "*11000000110000001*01001001001000001*10001100011000011",
/*  31*/  "*11111111001001011*00000110000111101*11111111001100000",
/*  32*/  "*00101011001111110*01010100001010010*01101011111111110",
/*  33*/  "*10011111000011100*00110111110000000*00101010010010010",
// /*  34*/  "*11101111110111111*11111110111011110*11001011100011010",
/*  35*/  "*10101010000111101*01010000000000010*00000001000000010",
/*  36*/  "*01101101101110000*11001011001010000*01000001011111111",
/*  37*/  "*01010010001001101*00000101011111110*01000010000111111",
/*  38*/  "*01100100001011100*01111000100100011*00010000100111001",
/*  39*/  "*01011000000111100*00011010001110100*10011111001000001",
/*  40*/  "*01010000111110111*00010001111100001*11100110011010010",
/*  41*/  "*00000111100011100*10101011010001100*10101111000110111",
/*  42*/  "*11111011111000000*01000011100000001*11111101100100010",
/*  43*/  "*01010101110100001*11011111111010110*01100001000110000",
/*  44*/  "*00010101010001110*00111010001111111*00001111101001111",
/*  45*/  "*11000101100000100*00001101011111110*01111011111101000",
/*  46*/  "*10010110000010000*00000011011101000*10010010111010011",
/*  47*/  "*10011011100111001*01100010000010000*00000011111010001",
/*  48*/  "*00111100011110010*11100000010100011*00111100000101100",
/*  49*/  "*11101110011011011*00100011111111001*11011001101111101",
/*  50*/  "*00011110100110001*01011100111101110*00000111011110111",
/*  51*/  "*10001011001010010*11111111100000011*11100111100110100",
/*  52*/  "*00000110101110001*00000110111010000*11111111111110101",
/*  53*/  "*01010111110100101*11111000101101000*01001000000100011",
/*  54*/  "*10110100101111011*11001111100010111*00011100011111000",
/*  55*/  "*01010010110000000*01110111001111111*11110100111001011",
/*  56*/  "*01100010011111111*10000111110101101*00010111111010110",
/*  57*/  "*10011101111111110*00000000101010100*11111010001111001",
/*  58*/  "*10011100001011100*11111001000010010*00001111110110011",
/*  59*/  "*01100011101000111*00010001100011111*01110010110110001",
/*  60*/  "*00111000010011111*11010111110000001*01111100110011110",
/*  61*/  "*10101010101111101*10010000000110000*11110001010001001",
/*  62*/  "*00010000011101110*01110000001010001*10010100110001111"};

int prop2(char *s)
// takes the concatenation of three strings and tests property 2 on it.
{
    int res[3*LEN+1]; 
   
    int l,r,k,i,j;

    for(l=0;l<LEN;l++)if(s[l]=='*'){

       EDstardiagprefnonstraight(s+l,3*LEN-l,s+l,3*LEN-l,res);

       int starcnt=0;

       for(r=l;r<3*LEN;r++)if(s[r]=='*'){

         starcnt++;
//         if(res[r-l+1] != EDstarnonstraight(s+l,r-l+1,s+l,r-l+1))printf("ED*NS=%d  %d    *cnt=%d\n",res[r-l+1], EDstarnonstraight(s+l,r-l+1,s+l,r-l+1),starcnt);
         if(res[r-l+1] < starcnt){

             printf("Problem Prop 2: [%d,%d] ED*NS=%d\n",l,r,res[r-l+1]);
             return -1;
         }
 

       }

   }

   return 1;

}


int testproperty2(int cnt)
{

    int i,j,k;

    char s[3*LEN+1];

    for(i=0;i<cnt;i++){
        for(j=0;j<cnt;j++)
            for(k=0;k<cnt;k++){

       strcpy(s,      colstr[i]);
       strcpy(s+LEN,  colstr[j]);
       strcpy(s+2*LEN,colstr[k]);
    
       if(prop2(s)==-1){

          printf("|  %3d %s\n",i,colstr[i]);
          printf("|  %3d %s\n",j,colstr[j]);
          printf("+- %3d %s\n",k,colstr[k]);fflush(stdout);
       }

    }
    printf("Prop 2: cycles done %d %d %d\n",i,j,k);fflush(stdout);
    }

    printf("Testing prop 2 done\n");
    return 0;
      
}

int prop4(char *s)
// takes the concatenation of three strings and tests property 4 on it.
{
    int res[3*LEN+1]; 
   
    int l,r,k,i,j;

    for(l=0;l<LEN;l++){

       EDstarprefnonstraight(s+l,3*LEN-l, s+l,2*LEN-l,res);


       for(r=2*LEN;r<3*LEN;r++){
//           printf("%d %d %d  %f \n",l,r,res[r-l+1],ALPHA * (2*LEN - l));

           if(res[r-l+1] < ALPHA * (2*LEN - l)){

                printf("Problem Prop 4: [%d,%d] vs [%d,%d] ED*NS=%d  thr=%f\n",l,2*LEN,l,r,res[r-l+1], ALPHA * (2*LEN - l));
                printf("      ");for(i=0;i<2*LEN-l;i++)printf("%c",s[l+i]);printf("\n");
                printf("  vs  ");for(i=0;i<r-l+1;i++)printf("%c",s[l+i]);printf("\n");
                return -1;
          }
       }
 
   }

   return 1;

}


int testproperty4(int cnt)
{

    int i,j,k;

    char s[3*LEN+1];
    char sr[3*LEN+1];

    for(i=0;i<cnt;i++){
        for(j=0;j<cnt;j++)
            for(k=0;k<cnt;k++){

       strcpy(s,      colstr[i]);
       strcpy(s+LEN,  colstr[j]);
       strcpy(s+2*LEN,colstr[k]);
    
       if(prop4(s)==-1){

          printf("|  %3d %s\n",i,colstr[i]);
          printf("|  %3d %s\n",j,colstr[j]);
          printf("+- %3d %s\n",k,colstr[k]);fflush(stdout);
       }

       reversestr(s,3*LEN,sr);

       if(prop4(s)==-1){
          printf("|  reversed\n");
          printf("|  %3d %s\n",i,colstr[i]);
          printf("|  %3d %s\n",j,colstr[j]);
          printf("+- %3d %s\n",k,colstr[k]);fflush(stdout);
       }
       

    }
    printf("Prop 4: cycles done %d %d %d\n",i,j,k);fflush(stdout);
    }

    printf("Testing prop 4 done\n");
    return 0;
      
}

int prop3b1(char *s)
// takes the concatenation of three strings and tests property 3b1 on it.
{
    int l,r,k,i,j;

    for(l=0;l<LEN;l++){

       int res=EDstarminrangenonstraight(s+l,LEN+LEN/4, s,LEN, l, LEN, -l);

           if(res < ALPHA*LEN){

                printf("Problem Prop 3b1: [%d,%d] vs [%d,%d] ED*NS=%d  thr=%f\n",l,LEN+l+res,1,LEN, res, ALPHA * LEN);
                printf("      ");for(i=0;i<LEN+res;i++)printf("%c",s[l+i]);printf("\n");
                printf("  vs  ");for(i=0;i<LEN;i++)printf("%c",s[i]);printf("\n");
                return -1;
          }
       
 
   }

   return 1;

}

int prop3b2(char *s)
// takes the concatenation of three strings and tests property 3b2 on it.
{
    int l,r,k,i,j;

    for(l=0;l<LEN;l++){

       int res=EDstarminrangenonstraight(s+l,LEN+LEN/4, s+LEN,LEN, 0, LEN, LEN-l);

           if(res < ALPHA*LEN){

                printf("Problem Prop 3b2: [%d,%d] vs [%d,%d] ED*NS=%d  thr=%f\n",l,LEN+l+res,1,LEN,res, ALPHA * LEN);
                printf("      ");for(i=0;i<LEN+res;i++)printf("%c",s[l+i]);printf("\n");
                printf("  vs  ");for(i=0;i<LEN;i++)printf("%c",s[LEN+i]);printf("\n");
                return -1;
          }
       
 
   }

   return 1;

}


int prop3b(char *s)
// takes the concatenation of three strings and tests property 3b on it.
{
   if(prop3b1(s)!=-1)return prop3b2(s);

   return -1;
}


int testproperty3b(int cnt)
{

    int i,j,k;

    char s[3*LEN+1];

    for(i=0;i<cnt;i++){
        for(j=0;j<cnt;j++)
            for(k=0;k<cnt;k++){

       if((i==j)||(j==k)||(i==k))continue;

       strcpy(s,      colstr[i]);
       strcpy(s+LEN,  colstr[j]);
       strcpy(s+2*LEN,colstr[k]);
    
       if(prop3b(s)==-1){

          printf("|  %3d %s\n",i,colstr[i]);
          printf("|  %3d %s\n",j,colstr[j]);
          printf("+- %3d %s\n",k,colstr[k]);fflush(stdout);
       }


    }
    printf("Prop 3b: cycles done %d %d %d\n",i,j,k);fflush(stdout);
    }

    printf("Testing prop 3b done\n");
    return 0;
      
}



int main(int argc, char* argv[])
{

   int i,j,k,l,r,offs;
   int testres;
   int rep=0;
   int remcnt;

   cnt=61;

   if(LEN!=strlen(colstr[0])){printf("LEN collection mismatch!!!\n");return 0;}

   printf("Testing collection LEN=%d RATE=%d ALPHA=%0.3f COLSIZE=%d\nCollection:\n", LEN, RATE, ALPHA, cnt);fflush(stdout);

   for(i=0;i<cnt;i++)printf("%3d %s\n",i, colstr[i]);

   printf("Testing properties:\n");

   testproperty2(cnt);
   testproperty4(cnt);
   testproperty3b(cnt);


   return 0;
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

int maincol(int argc, char* argv[])
{

   int i,j,k,l,r,offs;
   int testres;
   int rep=0;
   int remcnt;

   printf("LEN=%d RATE=%d RANKFRAC=%f\n", LEN, RATE, RANKFRAC);

   stringstat(10000);
   
   for(remcnt=offs=cnt=0;cnt<COLSIZE;){
      randombiasstarstring(collectionstr[cnt]);
      
      testres=1;
      for(j=0;(j<cnt)&&(testres>0);j++)testres=test(collectionstr[(j+offs)%cnt],collectionstr[cnt]);

      if(testres==1){
         printf("\n%4d %s\n",cnt,collectionstr[cnt]);fflush(stdout);
         cnt++;rep=0;collfail[cnt]=0;remcnt++;

         int maxfail=0;
         int maxk;
         for(k=0;k<cnt;k++){ printf("%ld ",collfail[k]); if(maxfail<collfail[k]){maxk=k; maxfail=collfail[k];} } printf("\n");
         if(remcnt>=10){
           printf("removing %d\n",maxk);
           strcpy(collectionstr[maxk],collectionstr[cnt-1]);
           collfail[maxk]=0;
           cnt--;remcnt=0;
           for(k=0;k<cnt;k++)printf("col[%4d]=\"%s\";\n",k,collectionstr[k]);fflush(stdout);
         }
      }else{
        offs=(1+j+offs)%cnt;
        collfail[offs]++;
        if((rep++%100)==0){printf("+");fflush(stdout);}
      }
  }

   return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////


int OKforcollection(int64 xint, int len)
{
    char x[MAXLEN], y[MAXLEN];
   
    inttostr(xint, x, len);

    if(slantedLCS(x,len,x,len)>0.85 * len)return 0;   // test self-similarity

    for(int i=0;collection[i]!=0;i++){
    
       inttostr(collection[i], y, len);
       
//       printf(". LCS=%d  x=%s   y=%s\n",LCS(x,len,y,len),x,y);
       
       if(LCS(x,len,y,len)>0.85 * len)return 0;       // test againsta all elements in the collection
    
    }
    
    return 1;

}

int mainx(int argc, char* argv[])
{

   char x[MAXLEN], y[MAXLEN];
   int64 xint;
   int n=1000,i,rept,ed;

#define REPET 3000

 
  for(i=0;i<10000;i++){
     collection[i]=0;

     do{
        xint=random64();
     }while(OKforcollection(xint,LEN)==0);
     
     collection[i]=xint;
     
     inttostr(xint,x,LEN);
     printf("%4d %s\n",i,x);fflush(stdout);
  }

   return 0;
}

int mainy(int argc, char* argv[])
// computing average LCS for random strings
{
   char x[MAXLEN], y[MAXLEN];
   int64 xint;
   int n=1000,i,j,rept,ed;

//   char *x1="abababab";
//   char *y1="bababab*";
//   printf("%d %s %s \n",EDstar(x1,strlen(x1),y1,strlen(y1)),x1,y1); return 0;

   srandom32(10); 

   for(n=49;n<50;n+=1){

      for(rept=ed=0;rept<REPET;rept++){
         randomstring(x,n);randomstring(y,2*n);
         strstar(x,n,RATE);strstar(y,2*n,RATE);
         // ed+=slantedLCS(x,strlen(x),x,strlen(x));
         ed+=EDstarmin(x,strlen(x),y,strlen(y));

         for(j=-10;j<10;j++)printf("%d %d   ",j,EDstar(x,strlen(x),y,strlen(x)+j));      
         printf("min %d\n",EDstarmin(x,strlen(x),y,strlen(y)));

         printf("%i %i %f %s %s\n",n,rept,((float)ed)/(n*(rept+1)),x,y);fflush(stdout);
      }
   }
   
  return 0;
}


int mainS(int argc, char* argv[])
// computing average LCS for random strings
{
    int i,j,k,l,r;

    char x[LEN+1],y[LEN+1];
    char xr[LEN+1],yr[LEN+1];

    while(j<1000000){

        randombiasstarstring(x);
        reversestr(x,LEN,xr);
        randombiasstarstring(y);
        reversestr(y,LEN,yr);

        if(EDstar(x,LEN,y,LEN)!=EDstar(x,LEN,y,LEN)){ printf("%d %d %d %d\n",EDstar(x,LEN,y,LEN),EDstar(x,LEN,y,LEN),EDstar(xr,LEN,yr,LEN),EDstar(xr,LEN,yr,LEN) );fflush(stdout);}

        if((++j%10000)==0)printf("%d\n",j);

   }

  return 0;
}
