# CCD
Cyclic Coordinate Decent ++ (CCD++) algorithm for matrix factorization, paralleled by MPI.

## 1 Description

### 1.1 Background

Matrix factorization is to decompose the sparse matrix R into two matrix U and V, and enables R≈UV. This problem can be solved by some machine learning method such as ALS, SGD, CCD, etc. CCD++ is a variant of CCD, proposed in [1]. CCD++ can be easily parallized and suitable for training large scale dataset. EsPa is a partition method originally used in CCD++, which may cause load imbalance and lowering the parallel efficiency. BaPa is an advanced partition method which effectively alleviates load imbalance problem and usually makes the program run faster than EsPa. To get more information about EsPa and BaPa, see our paper [2].

### 1.2 MPI implementation

Data transmission is implemented by blocking sending (MPI_Send) and non-blocking receiving (MPI_Irecv and MPI_Wait), which can effectively reduce the communication cost.

### 1.3 How to use our program

Before running our program, you need to prepare 2 datasets, whose data need to be arranged in specific format (see Section 2.1). Then you need to compile the program. Some marco values are required to be specified by adding them to the CFLAG argument, or you will receive compilation error (see Section 2.2). To run the program, you need to specify the number of MPI tasks before submit it to the cluster (see Section 2.3).

The program will generate 5 files. Here are their names and descriptions:

|NAME	|		DESCRIPTION |
| --- | --- |
|blockfile_row	|	Record the nonzeros contained in each 2D block of original dataset. |
|blockfile_col	|	Record the nonzeros contained in each 2D block of transposed dataset. |
|r_rtag.txt		|	Record the number of nonzeros of each 2D block of original dataset. |
|r_ctag.txt		|	Record the number of nonzeros of each 2D block of transposed dataset. |
|record.txt		|	Record the RMSE result of each iteration, "loop count" and "average elapsed time per iteration". |

"blockfile_row" and "blockfile_col" are the partition results of the datasets. "record.txt" records the performance results of our program.

The program currently does not output anything more than these files. If you want to see more details, such as the training result of U, V and R, you need to write print code yourself.


## 2 Getting Started

### 2.1 Prepare Datasets

You need to prepare 2 datasets as the input: the original dataset and its transposed dataset.

* **The original dataset**

The sparse matrix R should be stored in COO format, i.e, each nonzero is stored as a triple of (row index | column index | value). Details of the three elements are as follows:

|NAME			    | TYPE			|    RANGE       |
| --- | :---: | --- |
|row index	  | integer	  |		1 to USERS   |
|column index | integer	  |		1 to ITEMS   |
|value			  | double		| 	>=0          |

The nonzeros should be sorted by row index in ascending order. For those nonzeros having the same row index, they should be sorted by column index in ascending order.

The user ids of original dataset may be incontinuous. For example, they may be (1, 3, 10, 20...) We made a mapping of (1, 3, 10, 20...)->(1, 2, 3, 4...) so that the user ids could be continuous. Such mapping is also used to deal incoutinuous item ids.

* **The transposed dataset**

The transposed dataset is derived from the original dataset, by implementing the following 3 steps:
1. Swap the value of "row index" and "column index" for each nonzero. For example, the triple (1  6  5) in original dataset corresponds to the triple (6  1  5) in transposed dataset.
2. Sort the nonzeros by row index in ascending order.
3. Sort the nonzeros of the same row index by column index in ascending order.

We have provided an original dataset "movielens1m.dat" and its transposed dataset "movielens1m_trans.dat" for example.

### 2.2 Complie

Severel marcos are required to specify when compiling our program. Here are the lists of these marcos:

|  MACRO NAME        |   MACRO VALUE TYPE  |   DESCRIPTION |
| --- | :---: | --- |
|      ROW           |        integer      | The number of row partitioning. Matrix R is partitioned into ROW * COL blocks, where COL = ROW. |
|      USERS         |        integer      | The number of users. |
|      ITEMS         |        integer      | The number of items. |
| DATASET_ROW_PATH   |        string       | The path of your original dataset. |
| DATASET_COL_PATH   |        string       | The path of your transposed dataset. |
|ALGO_BAPA ALGO_ESPA |                     | The partition algorithm. Choose either BAPA or ESPA. |

When using "make" to compile the program, these marcos can be specified in the format 

```
CFLAGS+=-D<macro name(=macro value)>
```

Here is a compilation example:

```
make CFLAGS+=-DROW=12 CFLAGS+=-DUSERS=6040 CFLAGS+=-DITEMS=3952 CFLAGS+=-DALGO_BAPA CFLAGS+=-DDATASET_ROW_PATH=./movielens1m.dat CFLAGS+=-DDATASET_COL_PATH=./movielens1m_trans.dat
```

which can also be written as

```
make CFLAGS+="-DROW=12 -DUSERS=6040 -DITEMS=3952 -DALGO_BAPA -DDATASET_ROW_PATH=./movielens1m.dat -DDATASET_COL_PATH=./movielens1m_trans.dat"
```

### 2.3 Run

The program need an extra task to collect the results from other tasks. The extra task does not paticipate in parallel computation.  Remember that you have already compiled the program with a specific ROW (see Section 2.2). Now, if you want to run the program, you need to submit it with ROW+1 tasks.

The MPI program can be executed by specifying the task number when submitted to the cluster. On Tian-He2 platform, we use "yhrun" command to submit MPI program:

```
yhrun -N <minimum node number> -n <minimum task number> ./CCD.exe
```

For example, if you compile the program by setting CFLAGS+=-DROW=12, the submitting command becomes:

```
yhrun -N 13 -n 13 ./CCD.exe
```

It is OK if the "-N" argument less than 13, since a node may have several CPUs, and each CPU may have several cores. The minimum requirement is to ensure the task number not below the core number used, since one task is executed independently by one core.

Empirically, when "minimum node number" equals to "minimum task number", i.e., each node run only one task, the computation time is the lowest. If 2 or more tasks run on one node, they may scramble for computation resources and thus causing longer computation time.

If not on Tian-He2 platform, the command may be different.

## 3 Future Work
Here are some tips to further optimize our program, which have not been fully tested at present. We take them as our future work.

* Store data in lower precision
Our program stores the values in double precision. Actually high precision is not necessary in Machine Learing area. Nowadays many machine learning software incline to store data in single precision or half-single precision, in order to make the program run faster.

* Use optimizing argument 
Try "-O3" argument when compiling the program, e.g., adding "CFLAGS+=-O3" to the compilation example in Section 3.1. The argument "-O3" tells the compiler to produce a more optimized assemble code, which enables the program to run faster but probably with a slight precision loss.


## References:
[1] H.-F. Yu, C.-J. Hsieh, S. Si, and I. Dhillon, “Scalable coordinate descent approaches to parallel matrix factorization for recommender systems,” in 2012 IEEE 12th International Conference on Data Mining. IEEE, 2012, pp. 765–774.

[2] R. Guo et al., "BaPa: A Novel Approach of Improving Load Balance in Parallel Matrix Factorization for Recommender Systems," in IEEE Transactions on Computers, doi: 10.1109/TC.2020.2997051.
