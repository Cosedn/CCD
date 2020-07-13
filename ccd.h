#ifndef CCD_H
#define CCD_H

#ifndef ROW		// ROW is the row partitioning number, the matrix R is partitioned into
			// ROW*COL blocks where COL=ROW.
#error "please input the row partitioning number, e.g., CFLAGS+=-DROW=12\n"
#else
#define COL ROW
#endif

#define ATTR 5		// ATTR is the number of attributes of matrix U and V. The size of U is
			// USERS*ATTR, and the size of V is ITEMS*ATTR.
#define MAX_ITER 10000	// MAX_ITER is the maximum number of iterations. The training will
			// terminate when either of the following two conditions is met: (1) the
			// RMSE of current iteration is greater than the RMSE of previous iteration;
			// (2) the number of iterations exceeds MAX_ITER.
#define BUF 500
#define MAXIMUM 100000000.0

#ifdef MOVIELENS_100K
#define USERS 943	// USERS is the number of users.
#define ITEMS 1682	// ITEMS is the number of items.
#elif defined MOVIELENS_1M
#define USERS 6040
#define ITEMS 3952
#elif defined MOVIELENS_10M
#define USERS 69878
#define ITEMS 10677
#elif defined MOVIELENS_20M
#define USERS 138493
#define ITEMS 26744
#elif defined NETFLIX||NETFLIX2
#define USERS 480189
#define ITEMS 17770
#else
#ifndef USERS
#error "please input the number of users of your dataset, e.g., CFLAGS+=-DUSERS=6040\n"
#endif
#ifndef ITEMS
#error "please input the number of items of your dataset, e.g., CFLAGS+=-DITEMS=3952\n"
#endif
#endif

#ifndef DATASET_ROW_PATH	// DATASET_ROW_PATH is the path of your original dataset.
#error "please input the path of your original dataset, e.g., CFLAGS+=-DDATASET_ROW_PATH=./movielens1m.dat\n"
#endif

#ifndef DATASET_COL_PATH	// DATASET_COL_PATH is the path of the transposed dataset.
#error "please input the path of your transposed dataset, e.g., CFLAGS+=-DDATASET_COL_PATH=./movielens1m_trans.dat\n"
#endif

#define STR1(R) #R
#define STR2(R) STR1(R)  //double stringfication

#if defined ALGO_BAPA||ALGO_ESPA
#else
#error "Please select a partition algorithm. Input CFLAGS+=-DALGO_BAPA to use BAPA, or input CFLAGS+=-DALGO_ESPA to use ESPA\n"
#endif

typedef struct
{
	int user;
	int item;
	double rating;
}RatingMatrix;

typedef struct
{
	int row;
	int col;
	int entries;
}Block;

typedef struct
{
	int user;
	int item;
	int rating;
}RatingMatrixTemp;

int Input_File(char *, int *, int *);
void Make_Block(int *, int *, int *, int *, int, int, double);
void File_Block_Get(Block, char *, RatingMatrix *);
void File_Separate(char *, char *, int *, int *, Block *, Block *);
void Store_Rtag(char *, Block *);
void Get_Rtag(char *, Block *);

#endif
