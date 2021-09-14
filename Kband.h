typedef short int bool;

bool insideband(int i, int j, int band); 
int mapping(int x, int y, int n);
void init_allline(int n, int band, int *p);
double protein_score(char a, char b, double **mtx);
int query_place(int x, int y, int band, int *s);
int matrix_set(int x, int y, int band, double val, double *mtx, int *place_seq);
double matrix_query(int x, int y, int band, double *mtx, int *place_seq);
int matrix_update(int x, int y, int band, double val, double *mtx, int *place_seq);
int matrix_set_INT(int x, int y, int band, int val, int *mtx, int *place_seq);
int matrix_query_INT(int x, int y, int band, int *mtx, int *place_seq);
