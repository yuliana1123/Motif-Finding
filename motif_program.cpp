#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define N 1000
#define sn 50

double max(double * ,int );
int mo_len = 0 ;

int main()
{
	srand(time(NULL));
	char first[N] = {0} ,m_len[1][N]={0};
	int num ,n1 ,n2 ;
	int len = 0 ,time=0 ;
	char all_seq[sn][N] = {0} ;
	int seq_n = 0;
	int mo_len ;

	FILE *fp = fopen("Q1.txt", "r"); //change the file name
	
	if(fp == NULL)
	{
		printf("cannot open file\n");
		return 1;
	}
	else{
		printf("read successful\n");
	}
	
	//get motif length
	while(fgets(first,sizeof(first),fp) ){
		
		len = 0;
		
		mo_len = 15 ;
		
		//printf("motif width = %d\n",mo_len);
		
		char *find ;
		find = strstr(first,"seq");
		
		if(find != NULL){
			
			fgets(all_seq[seq_n],sizeof(all_seq[seq_n]),fp);
			seq_n++;
		}
		
	}
	
	fclose(fp);
	
	int i ,j ;
	int seq_len ;
	int beg_pos ;
	
	
	for( i = 1 ; i <= sn ; i++ ){
		printf("seq%d = %s\n",i ,all_seq[i-1] );
	}

	seq_len = strlen(all_seq[0]);
	
	//printf("seq_len = %d\n",seq_len);	
	
	beg_pos = (int)( (rand()/(RAND_MAX + 1.0)) * (seq_len - 2 - mo_len + 1 + 1.0 - 0) );
	
	//printf("beg_pos = %d\n",beg_pos );
	
	//printf("mo_len = %d\n",mo_len );
	//printf("sn = %d\n",sn );
	
	char **all_mot , **ex_fir_mot;
	void test_pro (char **a_sep , double a_pro[][4] ,int);
	double beg_fun_score[mo_len][4] = {0};
	
	all_mot = new char *[sn];
	for( i = 0 ; i < sn ; i++ )
		all_mot[i] = new char [mo_len];
	
	
	for( i = 0 ; i < sn ; i++ ){
		for( j = 0 ; j < mo_len ; j++ ){
			all_mot[i][j] = all_seq[i][beg_pos] ;
			beg_pos++;
		}
		beg_pos -= mo_len ;
	}
	
	//printf("beg_pos = %d\n",beg_pos );
	
	/*
	for( i = 0 ; i < sn ; i++ ){
		printf("all_mot[%d] = ",i );
		for( j = 0 ; j < mo_len ; j++ )
			printf("%c ",all_mot[i][j] );
		printf("\n");
	}
	*/	
	
	test_pro (all_mot ,beg_fun_score ,mo_len);
	
	/*
	printf("\n");
	for( i = 0 ; i < mo_len ; i++ ){
		printf("beg_fun_score[%d] = ",i );
		for( j = 0 ; j < 4 ; j++ )
			printf("%lf ",beg_fun_score[i][j] );
		printf("\n");	
	}
	*/
	
	/*
	for( i = 0 ; i < sn - 1 ; i++ )
		for( j = 0 ; j < mo_len ; j++ )
			printf("ex_fir_mot[%d][%d] = %c\n",i ,j ,ex_fir_mot[i][j] );
	*/	
	
	 //A,T,C,G
	int cou_time = seq_len - 2 - mo_len + 1;
	char motif[mo_len] = {0};
	int mo_pos[sn] = {0};
	//void cou_pro (char **sep , double pro[][4] ,int);
	void each_score (char * ,char * ,int *,int ,int ,int);
	void find_mose (double m_pro[][4] ,char * ,int );
	int error_time (char * ,char **all_motif ,int );
	
	//decide
	printf("counting...\n");
	int keep_end(double pro[][4] ,int );
	int prime2 = 0;
	do{
		int prime = 0;
		do{
			int exp_seq = 0 ;
			int now_seq[sn-1] = {0};
			beg_fun_score[mo_len][4] = {0};
			do{
				
				
				
				//cou_pro(ex_fir_mot ,beg_fun_score ,mo_len);
				each_score (all_seq[exp_seq] ,all_mot[exp_seq] ,mo_pos ,exp_seq ,cou_time ,mo_len);
					
				exp_seq++;
				
			}while(exp_seq < sn);
			
			beg_fun_score[mo_len][4] = {0};
			test_pro (all_mot ,beg_fun_score ,mo_len);
			prime = keep_end(beg_fun_score ,mo_len );
			
			//printf("\nprime = %d\n",prime );
		
		}while(prime != 1);
	
		find_mose(beg_fun_score ,motif ,mo_len );
		
		
		
		prime2 = error_time (motif ,all_mot ,mo_len );
	
	}while(prime2 != 1);
	
	printf("\nmotif = ");
		for( i = 0 ; i < mo_len ; i++ )
			printf("%c ",motif[i] );
		printf("\n");
	
	
	for( i = 0 ; i < sn ; i++ ){
		printf("\nseq %2d  #%3d ",i ,mo_pos[i] );
		for( j = 0 ; j < mo_len ; j++ )
			printf("%c ",all_mot[i][j] );
		printf("\n");
	}
	
	/*
	printf("\n");
	for( i = 0 ; i < mo_len ; i++ ){
		printf("\nbeg_fun_score[%d] = ",i );
		for( j = 0 ; j < 4 ; j++ )
			printf("%lf ",beg_fun_score[i][j] );
		printf("\n");	
	}
	*/
	
	system("pause");
	
	for( i = 0 ; i < sn ; i++ )
		delete []all_mot[i];
	delete []all_mot;
	
	
	
	
	
	
	
	
	

	return 0;
}

void each_score (char *one_seq ,char *allmot ,int *mo_pos ,int exp,int cou_time ,int mo){
	
	int i ,j ,r ;
	double sum = 1 ,m ,bg = 1;
	
	r = (int)( (rand()/(RAND_MAX + 1.0)) * (cou_time + 1.0 - 0) );
	
	mo_pos[exp] = r ;
	
	for(i = 0 ;i < mo ;i++ )
		allmot[i] = one_seq[r+i] ;
	
		
	
	
	
}

double max(double *m,int len){
	
	double max ;
	max = m[0];
	for(int i = 0;i < len;i++ ){
		if(m[i] > max)
			max = m[i];
			
	}
	
	return max ;
	
}

void test_pro (char **a_sep , double a_pro[][4] ,int mo){
	
	int i ,j ;
	double sum[mo] = {0};
	int cou_num[mo][4]= {0} ;
	
	for( j = 0 ; j < sn ; j++ ){
		for( i = 0 ; i < mo ; i++ )
		{
			
			if(a_sep[j][i] == 'A' || a_sep[j][i] == 'a' )
				cou_num[i][0] += 1;
			else if(a_sep[j][i] == 'T' || a_sep[j][i] == 't' )
				cou_num[i][1] += 1;
			else if(a_sep[j][i] == 'C' || a_sep[j][i] == 'c' )
				cou_num[i][2] += 1;
			else
				cou_num[i][3] += 1;		
			
		}
		
	}
	
	for( i = 0 ; i < mo ; i++ )
		for( j = 0 ; j < 4 ; j++ )
			a_pro[i][j] = 0.5 ;
	
	for( i = 0 ; i < mo ; i++ ){ 
		for( j = 0 ; j < 4 ; j++ ){
			
			if(j == 0)
				a_pro[i][j] += cou_num[i][j] ;
			else if(j == 1)
				a_pro[i][j] += cou_num[i][j] ;
			else if(j == 2)
				a_pro[i][j] += cou_num[i][j] ;
			else
				a_pro[i][j] += cou_num[i][j] ;
			
		}
	}

	for( i = 0 ; i < mo ; i++ ){ 
		for( j = 0 ; j < 4 ; j++ ){
			
			sum[i] += a_pro[i][j] ;	
	
		}
	}
	/*
	for( i = 0 ; i < mo ; i++ )
		printf("sum[%d] = %lf\n",i ,sum[i] );
	*/
	
	
	for( i = 0 ; i < mo ; i++ )
		for( j = 0 ; j < 4 ; j++ )
				a_pro[i][j] /= sum[i] ;
	/*
	printf("\n");
	for( i = 0 ; i < mo ; i++ ){
		printf("pro[%d] = ",i );
		for( j = 0 ; j < 4 ; j++ )
			printf("%lf ",a_pro[i][j] );
		printf("\n");	
	}		
	*/	
}


int keep_end(double pro[][4] ,int mo){
	
	int i ,j ,prime = 0 ;
	int count0 = 0 ,count1 = 0 ,count2 = 0 ,count3 = 0 ;
	
	for(i = 0;i < mo;i++ ){
		for(j = 0;j < 4;j++ ){
			
			if(pro[i][j] >= 0.7)
				count0++;
			
			if(pro[i][j] >= 0.6)
				count1++;
				
			if(pro[i][j] >= 0.5)
				count2++;
				
			if(pro[i][j] >= 0.4)
				count3++;
		}
	}
	//if(count0 >= mo / 4.5 && count1 >= mo / 1.8 && count2 >= mo / 1.29 )
	//if(count0 >= mo / 4.5 && count1 >= mo / 2.25 && count2 >= mo / 1.29 && count3 >= mo / 1.25 )
	if(count3 >= mo / 2.75 )
		return 1;
	else
		return 0;
	
	
}

void find_mose (double m_pro[][4] ,char *motif_seq ,int mo){
	
	double num = 0 ,b = 0;
	char base;
	int i ,j ;
	
	for(i = 0;i < mo;i++ ){
		b = 0;
		num = 0;
		for(j = 0;j < 4;j++ ){
			
			if(m_pro[i][j] > num){
				b = j;
				num = m_pro[i][j];
			}
				
			
		}
		//printf("b = %lf\n",b);
		if (b == 0)
			base = 'A';
		else if (b == 1)
			base = 'T';
		else if (b == 2)
			base = 'C';
		else
			base = 'G';
			
		motif_seq[i] = base ;
		
	}
	
	
	
}

int error_time (char *motif_seq ,char **all_motif ,int mo){
	
	int i ,j ;
	double count = 0;
	double max_err = 0.8 * mo ;
	
	for(i = 0;i < sn;i++ ){
		
		count = 0;
		
		for(j = 0;j < mo;j++ ){
			
			if(all_motif[i][j] != motif_seq[j])
				count++;
			
		}
		if(count > max_err)
			return 0 ;
	}
	
	return 1 ;
	
}
