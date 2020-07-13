#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "ccd.h"

typedef struct _Head
{
	struct _Node *node;
	struct _Node *current;
}Head;

typedef struct _Node
{
	RatingMatrix buffer[BUF];
	int size;
	struct _Node *next;
}Node;

void Init_List(Head *head)
{
	Node *node;
	node = (Node *)calloc(1, sizeof(Node));
	node->size = 0;
	node->next = NULL;
	head->node = node;
	head->current = node;
}

void Add_To_List(Head *head, RatingMatrix *rm)
{
	if(head->current->size >= BUF - 1)
	{
		Node *node = NULL;
		if((node = (Node *)calloc(1, sizeof(Node))) == NULL)
		{
			printf("memory allocation failed!\n");
			exit(0);
		}
		node->size = 0;
		node->next = NULL;
		head->current->next = node;
		head->current = head->current->next;

		head->current->buffer[head->current->size] = *rm;
		head->current->size++;
	}
	else
	{
		head->current->buffer[head->current->size] = *rm;
		head->current->size++;
	}
}

void Destroy_List(Head *head)
{
	Node *temp;
	head->current = head->node;
	while(head->current != NULL)
	{
		temp = head->current->next;
		free(head->current);
		head->current = temp;
	}
}

/*void itoa(int num, char* str)
{
        int i=0,j;
        int t;
        char str1[10];
        while(1)
        {
                if(num==0) break;
                t=num-num/10*10;
                num=(num-t)/10;
                str1[i]=(char)(t+48);
                i++;
        }
        for(j=0;j<i;j++) str[j]=str1[i-j-1];
        str[j]='\0';
}*/

int Input_File(char* str, int *u_entries, int *v_entries)
{
	FILE *fp;
	int user,item;
	double rating;
	int num=0;
	if((fp=fopen(str,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
    	while(fscanf(fp,"%d %d %lf",&user,&item,&rating)!=EOF)
	{
		u_entries[user-1]++;
		v_entries[item-1]++;
		num++;
	}
	fclose(fp);
	return num;
}

#ifdef ALGO_BAPA
void Make_Block(int *u_entries, int *v_entries, int *u_block, int *v_block, int user, int item, double aver)
{
/*	int i;
	double row_number,col_number;
	row_number=(double)user/ROW;
	col_number=(double)item/COL;
	u_block[0]=0;
	v_block[0]=0;
	u_block[ROW]=user;
	v_block[COL]=item;
	for(i=1;i<ROW;i++) u_block[i]=(int)(i*row_number);
	for(i=1;i<COL;i++) v_block[i]=(int)(i*col_number);*/
	int i, j;
	double t0, t1, sum, aver_adder;
	u_block[0] = 0;
	i = 1;
	j = 0;
	sum = 0;
	t0 = MAXIMUM;
	aver_adder = aver;
	while(j < user)
	{
		sum += u_entries[j];
		t1 = (sum > aver_adder)?(sum - aver_adder):(aver_adder - sum);
		if(t1 <= t0)
		{
			t0 = t1;
			j++;
		}
		else
		{
			u_block[i] = j;
			sum -= u_entries[j];
			i++;
			t0 = MAXIMUM;
			aver_adder += aver;
		}
		if(i == ROW)
		{
			u_block[i] = user;
			break;
		}
	}
	v_block[0] = 0;
	i = 1;
	j = 0;
	sum = 0;
	t0 = MAXIMUM;
	aver_adder = aver;
	while(j < item)
	{
		sum += v_entries[j];
		t1 = (sum > aver_adder)?(sum - aver_adder):(aver_adder - sum);
		if(t1 <= t0)
		{
			t0 = t1;
			j++;
		}
		else
		{
			v_block[i] = j;
			sum -= v_entries[j];
			i++;
			t0 = MAXIMUM;
			aver_adder += aver;
		}
		if(i == COL)
		{
			v_block[i] = item;
			break;
		}
		//printf("%d %d %d %.2f %.2f\n", j, v_entries[j], i, sum, aver_adder);
	}
}
#endif

#ifdef ALGO_ESPA
void Make_Block(int *u_entries, int *v_entries, int *u_block, int *v_block, int user, int item, double aver)
{
	int i;
	double row_number,col_number;
	row_number=(double)user/ROW;
	col_number=(double)item/COL;
	u_block[0]=0;
	v_block[0]=0;
	u_block[ROW]=user;
	v_block[COL]=item;
	for(i=1;i<ROW;i++) u_block[i]=(int)(i*row_number);
	for(i=1;i<COL;i++) v_block[i]=(int)(i*col_number);
}
#endif
	
void File_Block_Get(Block rl_tag, char *type, RatingMatrix *rm_buffer) //Get the corresponding blocks in files for each thread
{
	FILE *fp;
	int i;
	int sz = 0;
	int size = 0;
	int match;
	int num_block = 0;
	RatingMatrix temp;

	if(strcmp(type,"row")==0)
	{
		if((fp=fopen("blockfile_row","r"))==NULL)
                {
                        printf("cannot open this file!\n");
                        exit(0);
                }
		fscanf(fp, "%d", &match);
		if(match != ROW)
		{
			printf("This file does not match the program!\n");
			exit(0);
		}

		while(1)
		{
			if(fscanf(fp, "%d", &size) == EOF) break;
			if(rl_tag.row == num_block)
			{
				for(i = 0; i < size; i++)
				{
					fscanf(fp, "%d %d %lf", &rm_buffer[sz].user, &rm_buffer[sz].item, &rm_buffer[sz].rating);
					sz++;
				}
				break;
			}
			else
			{
				for(i = 0; i < size; i++)
					fscanf(fp, "%d %d %lf", &temp.user, &temp.item, &temp.rating);
				num_block++;
			}
		}
                fclose(fp);
	}
	else if(strcmp(type,"col")==0)
	{
                if((fp=fopen("blockfile_col","r"))==NULL)
                {
                        printf("cannot open this file!\n");
                        exit(0);
                }
		fscanf(fp, "%d", &match);
		if(match != COL)
		{
			printf("This file does not match the program!\n");
			exit(0);
		}

		while(1)
		{
			if(fscanf(fp, "%d", &size) == EOF) break;
			if(rl_tag.col == num_block)
			{
				for(i = 0; i < size; i++)
				{
					fscanf(fp, "%d %d %lf", &rm_buffer[sz].user, &rm_buffer[sz].item, &rm_buffer[sz].rating);
					sz++;
				}
				break;
			}
			else
			{
				for(i = 0; i < size; i++)
					fscanf(fp, "%d %d %lf", &temp.user, &temp.item, &temp.rating);
				num_block++;
			}
		}
                fclose(fp);	
	}
}

void File_Separate(char *str1, char *str2, int *u_block, int *v_block, Block *r_rtag, Block *r_ctag) //Separate the Rating Matrix into blocks and store them in files
{
	FILE *fp, *fq;
	int i, k;
	int r, c, p;
	int user, item;
	double rating;
	int size;
	RatingMatrix rm;

	Head head[ROW];
	int num[ROW];

/*---------------------row (original) file separate---------------------------*/

	if((fp=fopen(str1,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	if((fq=fopen("blockfile_row","w"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	size = ROW;

	for(i = 0; i < ROW; i++) Init_List(&head[i]);

	while(fscanf(fp,"%d %d %lf",&user,&item,&rating) != EOF)
	{
		r = user - 1;
		c = item - 1;
		for(i = 0; i <= ROW; i++) if(r < u_block[i]) break;
		p = i - 1;
		rm.user = r;
		rm.item = c;
		rm.rating = rating;
		Add_To_List(&head[i-1], &rm);
		r_rtag[p].entries++;
	}

	for(i = 0; i < ROW; i++)
	{
		num[i] = 0;
		head[i].current = head[i].node;
		while(head[i].current != NULL)
		{
			num[i] += head[i].current->size;
			head[i].current = head[i].current->next;
		}
	}

	fprintf(fq, "%d\n", size);
	for(i = 0; i < ROW; i++)
	{
		fprintf(fq, "%d\n", num[i]);
		head[i].current = head[i].node;
		while(head[i].current != NULL)
		{
			for(k = 0; k < head[i].current->size; k++)
				fprintf(fq,"%d    %d    %.2f\n",
					head[i].current->buffer[k].user,
					head[i].current->buffer[k].item,
					head[i].current->buffer[k].rating);
			head[i].current = head[i].current->next;
		}
	}

	Store_Rtag("r_rtag.txt",r_rtag);
	fclose(fp);
	fclose(fq);
	for(i = 0; i < ROW; i++) Destroy_List(&head[i]);

/*---------------------column (transposed) file separate-----------------------*/

	if((fp=fopen(str2,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	if((fq=fopen("blockfile_col","w"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	size = COL;

	for(i = 0; i < ROW; i++) Init_List(&head[i]);

	while(fscanf(fp,"%d %d %lf",&user,&item,&rating) != EOF)
	{
		r = user - 1;
		c = item - 1;
		for(i = 0; i <= COL; i++) if(r < v_block[i]) break;
		p = i - 1;
		rm.user = r;
		rm.item = c;
		rm.rating = rating;
		Add_To_List(&head[i-1], &rm);
		r_ctag[p].entries++;
	}

	for(i = 0; i < COL; i++)
	{
		num[i] = 0;
		head[i].current = head[i].node;
		while(head[i].current != NULL)
		{
			num[i] += head[i].current->size;
			head[i].current = head[i].current->next;
		}
	}

	fprintf(fq, "%d\n", size);
	for(i = 0; i < COL; i++)
	{
		fprintf(fq, "%d\n", num[i]);
		head[i].current = head[i].node;
		while(head[i].current != NULL)
		{
			for(k = 0; k < head[i].current->size; k++)
				fprintf(fq,"%d    %d    %.2f\n",
					head[i].current->buffer[k].user,
					head[i].current->buffer[k].item,
					head[i].current->buffer[k].rating);
			head[i].current = head[i].current->next;
		}
	}

	Store_Rtag("r_ctag.txt",r_ctag);
	fclose(fp);
	fclose(fq);
	for(i = 0; i < COL; i++) Destroy_List(&head[i]);
}

void Store_Rtag(char *str, Block *r_tag)
{
	FILE *fp=NULL;
	int size=ROW;
	int i;
	if((fp=fopen(str,"w"))==NULL)
	{
		printf("cannot write to this file!\n");
		exit(0);
	}
	for(i=0;i<size;i++) fprintf(fp,"%d %d %d\n",r_tag[i].row,r_tag[i].col,r_tag[i].entries);
	fclose(fp);
}

void Get_Rtag(char *str, Block *r_tag)
{
	FILE *fp=NULL;
	int size=ROW;
	int i;
	if((fp=fopen(str,"r"))!=NULL)
	{
		for(i=0;i<size;i++) fscanf(fp,"%d %d %d",&r_tag[i].row,&r_tag[i].col,&r_tag[i].entries);
	}
	fclose(fp);
}

