#include "stdio.h"
#include "stdlib.h"
#include <string.h> 
#include <math.h> 
#include "time.h"
#include <pthread.h>

#define MAX_GENE 22270

float quickES(short s1[], short s2[],int len,int sig);

//short Source[MAX_SAMPLE][MAX_GENE];  //source gene set
//short Mat[MAX_SAMPLE][MAX_GENE];     //gene lib
//float ES[MAX_SAMPLE][MAX_SAMPLE];    //result ES matrix 50G

short **Source;  //source gene set
short **Mat;     //gene lib
float **ES;    //result ES matrix 50G

struct argument
{
	int sig;   //signaturn len 
    int genelen;   //gene len	
	int sourcelen;  //source len 
	int matdown;  //gene lib index down limit
	int matup;  //gene lib index up limit
};

void *thread_fuc(void *arg)
{
    int i,j;
    struct argument *arg_thread = (struct argument *) arg;
	
	for (i=0; i<arg_thread->sourcelen; i++)
	{
		for(j=arg_thread->matdown; j<arg_thread->matup; j++)
		{
			ES[i][j] = quickES(Source[i],Mat[j],arg_thread->genelen,arg_thread->sig);
		}
	}
	
	printf("sig:%d\tglen:%d\tslen:%d\tmdown:%d\tmup:%d\n",arg_thread->sig,arg_thread->genelen,arg_thread->sourcelen,arg_thread->matdown,arg_thread->matup);
}

void thread_create(pthread_t thread[],struct argument arg[],int cores)
{
    int i,temp;
    //create thread
	for(i=0; i<cores; i++)
	{
		if((temp = pthread_create(&thread[i],NULL,thread_fuc,(void *)&arg[i]))!=0)
			printf("thread%d create fail!\n",i+1);
		else
			printf("thread%d create!\n",i+1);
	}
}

void thread_wait(pthread_t thread[],int cores)
{
    //wait threads over
	int i;
	for(i=0; i<cores; i++)
	{
		if(thread[i]!=0)
		{
			pthread_join(thread[i],NULL);
			printf("thread%d over!\n",i+1);
		}
	}
    
}

float quickES(short s1[], short s2[],int len,int sig)
{
	int isgsUp[MAX_GENE] ;  
	int isgsDown[MAX_GENE] ;

	float ScorehitUp[MAX_GENE] ;
	float ScoremissUp[MAX_GENE]  ;
	float ScorehitDown[MAX_GENE] ;
	float ScoremissDown[MAX_GENE] ;

	int indexS[MAX_GENE] ; //gene local index
	
	int flagUp=0,flagDown=0;
	int hitsumUp=sig,hitsumDown=sig;
	float maxUp,maxDown;
	int indexUp,indexDown;
	float ES[2];
        
	int i,j;
	
	for( i=0; i<2; i++){		
		
		memset(isgsUp, 0, len * sizeof(int) );
		memset(isgsDown, 0, len * sizeof(int) );  //Init Up¡¢Down

		
		if(i==0){ 
			//build index for s2
			for( j=0 ; j<len; j++)
				indexS[s2[j]] = j;
			//compute isgsUp ¡¢ isgsDown
			for( j=0 ; j<sig; j++)
			{
				isgsUp[indexS[s1[j]]] = 1;
				isgsDown[indexS[s1[len-j-1]]] = 1;
			}				
		}else{
			//build index for s1
			for( j=0 ; j<len; j++)
				indexS[s1[j]] = j;
			//compute isgsUp ¡¢ isgsDown
			for( j=0 ; j<sig; j++)
			{
				isgsUp[indexS[s2[j]]] = 1;
				isgsDown[indexS[s2[len-j-1]]] = 1;
			}				
		}

		ScorehitUp[0]=isgsUp[0];
		ScoremissUp[0]=1-isgsUp[0];
		ScorehitDown[0]=isgsDown[0];
		ScoremissDown[0]=1-isgsDown[0];
		maxUp=fabs(ScorehitUp[0]/hitsumUp-ScoremissUp[0]/(len-hitsumUp));
		indexUp=0;
		maxDown=fabs(ScorehitDown[0]/hitsumDown-ScoremissDown[0]/(len-hitsumDown));
		indexDown=0;
			
		for( j=1;j<len;j++){
			
			//compute up gene
			ScorehitUp[j] = isgsUp[j]+ScorehitUp[j-1];
			ScoremissUp[j] = (1-isgsUp[j])+ScoremissUp[j-1]; //cumsum
			if(fabs(ScorehitUp[j]/hitsumUp-ScoremissUp[j]/(len-hitsumUp))>maxUp){
				maxUp	=	fabs(ScorehitUp[j]/hitsumUp-ScoremissUp[j]/(len-hitsumUp));
				indexUp	=	j;
			}
				
			//compute down	gene
			ScorehitDown[j] = isgsDown[j]+ScorehitDown[j-1];
			ScoremissDown[j] = (1-isgsDown[j])+ScoremissDown[j-1];
			if(fabs(ScorehitDown[j]/hitsumDown-ScoremissDown[j]/(len-hitsumDown))>maxDown){
				maxDown	= fabs(ScorehitDown[j]/hitsumDown-ScoremissDown[j]/(len-hitsumDown));
				indexDown	=	j;
			}
		}

		float ESUp = ScorehitUp[indexUp]/hitsumUp-ScoremissUp[indexUp]/(len-hitsumUp);
		float ESDown = ScorehitDown[indexDown]/hitsumDown-ScoremissDown[indexDown]/(len-hitsumDown);
			
		ES[i] =  ( ESUp - ESDown )/2;
	}

	return (ES[0]+ES[1])/2;
}

int ReadLine(char path1[],char path2[],int* sourcelen,int *matlen,int *genelen)
{
	FILE *fp; 
	char StrLine[135168]; 
	int strlen;           
	char c[] = " ";
	int line,col;
	
	int i,j;
	
	//read Source
	if((fp = fopen(path1,"r")) == NULL) 
	{ 
		printf("filel error!\n"); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	*sourcelen = atoi(strtok(StrLine,c));
	*genelen = atoi(strtok(NULL,c));
	
	if(*genelen>1000)
		strlen = 135168; //level3
	else
		strlen = 6144; //L1000
	
	Source = (short **)malloc(sizeof(short *)*(*sourcelen));  //dynamical init Source gene set 
	for(i=0; i<(*sourcelen); i++)
		Source[i] = (short *)malloc(sizeof(short)*(*genelen));
	//malloc will be free automatic till program over £¬not function return
	
	line = 0;
	while (!feof(fp)) 
    { 
		fgets(StrLine,strlen,fp);  //read one line
		col = 0;
		if(line<(*sourcelen)&&col<(*genelen))   
			Source[line][col++] = atoi(strtok(StrLine,c));
		char *p = strtok(NULL,c);
		while(p)
		{
			if(line<(*sourcelen)&&col<(*genelen))
				Source[line][col++] = atoi(p);
			p = strtok(NULL,c); 
		}
		line++;
	} 
	
	
	fclose(fp);
	
	
	//read Mat
	if((fp = fopen(path2,"r")) == NULL) 
	{ 
		printf("file2 error!\n"); 
		return -1; 
	} 
	fgets(StrLine,256,fp);
	*matlen = atoi(strtok(StrLine,c));
	
	Mat = (short **)malloc(sizeof(short *)*(*matlen));  //dynamical init Mat gene lib 
	for(i=0; i<(*matlen); i++)
		Mat[i] = (short *)malloc(sizeof(short)*(*genelen));
	
	ES = (float **)malloc(sizeof(float *)*(*sourcelen));  //dynamical init ES Matrix 
	for(i=0; i<(*sourcelen); i++)
		ES[i] = (float *)malloc(sizeof(float)*(*matlen));
		
	line = 0;
	while (!feof(fp)) 
    { 
		fgets(StrLine,strlen ,fp); 
		col = 0;
		if(line<(*matlen)&&col<(*genelen))
			Mat[line][col++] = atoi(strtok(StrLine,c));
		char *p = strtok(NULL,c);
		while(p)
		{
			if(line<(*matlen)&&col<(*genelen))
				Mat[line][col++] = atoi(p);
			p = strtok(NULL,c); 
		}
		line++;
	} 
	
	fclose(fp);                
	
	return 0; 
}

void WriteResult(int sourcelen ,int matlen ,int path[])
{
	int i;	
	
	FILE *fp;
	if((fp=fopen(path,"wb"))==NULL)
	{
		printf("can not open the file to write!");
		exit(1);
	}
	
	fwrite(&sourcelen,sizeof(int),1,fp);
	fwrite(&matlen,sizeof(int),1,fp);
	for( i=0; i < sourcelen; i++)
		fwrite(ES[i],sizeof(float),matlen,fp);
	
	fclose(fp);
}

void main( int argc, char *argv[] )
{	

	int i,j;
	int corenum = atoi(argv[1]);
	int sig = atoi(argv[2]);
	
	int sourcelen,matlen,genelen; 
	
	
	clock_t start,finish;
	double duration;
	start=clock();
	
	ReadLine(argv[3],argv[4],&sourcelen,&matlen,&genelen);  //get para from file
	
	pthread_t *thread = (pthread_t *)malloc(sizeof(pthread_t)*corenum);
    struct argument *arg = (struct argument *)malloc(sizeof(struct argument)*corenum);
		
	//split mat
	if(matlen%corenum==0)
	{  
		int per = matlen/corenum;
		for( i=0; i<corenum; i++)
		{
			arg[i].sig = sig;
			arg[i].genelen = genelen;
			arg[i].sourcelen = sourcelen;
			arg[i].matdown = i*per;
			arg[i].matup = (i+1)*per;
		}		
	}
	else
	{
		int per = matlen/corenum+1;
		for( i=0; i<corenum; i++)
		{
			arg[i].sig = sig;
			arg[i].genelen = genelen;
			arg[i].sourcelen = sourcelen;
			arg[i].matdown = i*per;
			
			if(i<(corenum-1))
			{				
				arg[i].matup = (i+1)*per;
			}
			else
			{  
				arg[i].matup = matlen;
			}				 
		}
	}		
	
	finish=clock();
	duration=(double)(finish-start)/CLOCKS_PER_SEC;
	printf("pre-operate: %.4f s\n",duration); 
	
	
	printf("I'm main func,I'm creating threads!\n");
	thread_create(thread,arg,corenum);
    printf("I'm main func,I'm waiting threads over!\n");
    thread_wait(thread,corenum);
	

	start=clock();

	WriteResult( sourcelen , matlen);
	
	finish=clock();
	duration=(double)(finish-start)/CLOCKS_PER_SEC;
	printf("write file operate: %.4f s\n",duration); 

/*
	for(i=0; i<sourcelen; i++)
	{
		for(j=0; j<matlen; j++){
			printf("%f\t",ES[i][j]);
		}
		printf("\n");
	}
	
*/

}

