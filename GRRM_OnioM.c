//
// GRRM_OnioM.c:
// Copyright (C) 2020-2024 Kimichi Suzuki
// Author: Kimichi Suzuki <ki_suzuki@eis.hokudai.ac.jp>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Last modified: 2022/10/22 
//
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <sys/types.h>
# include <ctype.h>
# include <unistd.h>
//
#define BUFSIZE  500
// 
typedef struct MsMlib{
	int            njobs, natom, matom;
	char          *inp, *MsMqM, *MsMmM;
//
	double    En[4];
	double  **XYZ;      
	char    **HML;      
	char    **Typ;      
	char    **Lnk;      
	char    **Mod;      
	int      *Con;      
	int      *Opt;      
	double  **MSQ;      
	double  **Grd0;     
	double  **Grd1;     
	double  **Grd2;     
	double  **Grd3;     
//
} MsMlib;
//

int                         **int_Matrix(const int size_a, int size_b);
void                          free_int_Matrix(int **matrix);
double                      **double_Matrix(const int size_a, int size_b);
void                          free_double_Matrix(double **matrix);
double                       *double_Vector(const int size_a);
void                          free_double_Vector(double  *matrix);
void                         *char_Matrix(int size_a, int size_b);
void                          free_char_Matrix(char **matrix);
int                          *int_Vector(const int n);
void                          free_int_Vector(int *v);
double                     ***double_3D_malloc(int size_a, int size_b, int size_c);
void                          double_3D_dealloc(double ***matrix);

//
static void    get_OnioM_prm     (MsMlib *pt);
// 
static int     MsMonioM_ReL      (MsMlib *pt);
static int     exe_MsMonioM_ReL  (MsMlib *pt);
static int     get_MsMonioM_ReL  (MsMlib *pt);
//
static int     MsMonioM_MoL      (MsMlib *pt);
static int     exe_MsMonioM_MoL  (MsMlib *pt);
static int     get_MsMonioM_MoL  (MsMlib *pt);
// 
static int     MsMonioM_MoH      (MsMlib *pt);
static int     exe_MsMonioM_MoH  (MsMlib *pt);
static int     get_MsMonioM_MoH  (MsMlib *pt);
//
static int     MsMonioM_job      (MsMlib *pt);
static int     MsMonioM_out      (MsMlib *pt);
//
static int     engrad_MsMonioM   (MsMlib *pt);
static int     prep_Hessian      (MsMlib *pt);
static int     hess_OnioM        (MsMlib *pt);
static int     hessian_max       (int a, int b);
static int     hessian_min       (int a, int b);
static int     hessian_index     (int a, int b);

double          scale_factor     ( char Sq, char Sm );
double          get_factor       ( char *s_qm, char *s_link, char *s_mm);
static int      prep_model_str   ( MsMlib *pt);

static int     ok, check;

// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
int           main (int argc, char *argv[]){

		char  *a;
	//
	//if ( argc == 2 ){
	if ( argc == 4 ){

		struct MsMlib  pp;
		struct MsMlib *pt = &pp;

		pt->inp   = &argv[1][0];
		pt->MsMqM = &argv[2][0];
		pt->MsMmM = &argv[3][0];
		check     =  0;

			check = MsMonioM_job(pt);
		pt->HML       = char_Matrix  (pt->natom, 7);
		pt->Typ       = char_Matrix  (pt->natom,50);
		pt->Mod       = char_Matrix  (pt->natom,50);
		pt->Lnk       = char_Matrix  (pt->natom, 3);
		pt->XYZ       = double_Matrix(pt->natom, 3);
		pt->MSQ       = double_Matrix(pt->natom, 3);
		pt->Grd0      = double_Matrix(pt->natom, 3);
		pt->Grd3      = double_Matrix(pt->natom, 3);
		pt->Grd2      = double_Matrix(pt->natom, 3);
		pt->Grd1      = double_Matrix(pt->natom, 3);
		pt->Con       = int_Vector   (pt->natom);
		pt->Opt       = int_Vector   (pt->natom);

		get_OnioM_prm(pt);
		if         ( check == 0 && pt->njobs != 0 ) {
			check =  MsMonioM_ReL        (pt);
			check =  exe_MsMonioM_ReL    (pt);
			check =  get_MsMonioM_ReL    (pt);
			//
			check =  prep_model_str      (pt);
			//
			check =  MsMonioM_MoL        (pt);
			check =  exe_MsMonioM_MoL    (pt);
			check =  get_MsMonioM_MoL    (pt);
			//
			check =  MsMonioM_MoH        (pt);
			check =  exe_MsMonioM_MoH    (pt);
			check =  get_MsMonioM_MoH    (pt);
			//
			check =  engrad_MsMonioM     (pt);
			if ( pt->njobs == 2 ){
				check = prep_Hessian (pt);
				if ( check == 0 ) check = hess_OnioM   (pt);
			}
			check = MsMonioM_out(pt);
		}

		free_int_Vector   (pt->Opt);
		free_int_Vector   (pt->Con);
		free_double_Matrix(pt->MSQ);
		free_double_Matrix(pt->XYZ);
		free_double_Matrix(pt->Grd0);
		free_double_Matrix(pt->Grd1);
		free_double_Matrix(pt->Grd2);
		free_double_Matrix(pt->Grd3);
		free_char_Matrix  (pt->Lnk);
		free_char_Matrix  (pt->Mod);
		free_char_Matrix  (pt->Typ);
		free_char_Matrix  (pt->HML);


	}
	//
	return 0;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
static void     get_OnioM_prm ( MsMlib *pt){

	FILE   *fp;
	char   *out, s[BUFSIZE];
	int     len, i, a, b;
	double  x, y, z;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 11 + 1 + 3 + 1)*sizeof(char));
	sprintf(out, "%s_GauJOB.com.prm\0", pt->inp);

	fp   = fopen (out, "r");
	if  (fp){

		for  ( i=0; i<pt->natom; ++i ){
			pt->XYZ[i][0]  = 0.0;
			pt->XYZ[i][1]  = 0.0;
			pt->XYZ[i][2]  = 0.0;
			pt->Opt[i]     =   0;
			pt->Con[i]     =   0;
			fgets (s, BUFSIZE, fp);
			sscanf(s, "%s%d%lf%lf%lf%s%d", pt->Typ[i], &a, &x, &y, &z, pt->HML[i], &b);
			pt->XYZ[i][0]  = x;
			pt->XYZ[i][1]  = y;
			pt->XYZ[i][2]  = z;
			pt->Opt[i]     = a;
			pt->Con[i]     = b;
		}
		fclose(fp);
	}
	free(out);
			//
			pt->     En[0] = 0.0;
			pt->     En[1] = 0.0;
			pt->     En[2] = 0.0;
			pt->     En[3] = 0.0;
		for  ( i=0; i<pt->natom; ++i ){
			pt->Grd0[i][0] = 0.0; pt->Grd0[i][1] = 0.0; pt->Grd0[i][2] = 0.0;
			pt->Grd1[i][0] = 0.0; pt->Grd1[i][1] = 0.0; pt->Grd1[i][2] = 0.0;
			pt->Grd2[i][0] = 0.0; pt->Grd2[i][1] = 0.0; pt->Grd2[i][2] = 0.0;
			pt->Grd3[i][0] = 0.0; pt->Grd3[i][1] = 0.0; pt->Grd3[i][2] = 0.0;
		}
		//
		for ( i = 0; i < pt->natom; ++i){
			if      ( strstr(pt->HML[i], "H")  ||  ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ) ++pt->matom;
		}
		//
}
// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
static int    MsMonioM_job (MsMlib *pt){

	FILE    *fp, *bp;
	char     sp[BUFSIZE], sq[BUFSIZE], *inpA, *inpB;
	double   rx, ry, rz;
	int      len;
	int      i, j, ij;
	char       p1[3], p6[3], p7[3];
	double     p3, p4, p5;
	int        p2, p8;

	//
	len  = strlen(pt->inp);
	inpA = (char*)malloc((len + 11 + 1)*sizeof(char));
	sprintf(inpA, "%s_GauJOB.com\0", pt->inp);
	//
	inpB = (char*)malloc((len + 11 + 4 + 1)*sizeof(char));
	sprintf(inpB, "%s_GauJOB.com.prm\0", pt->inp);
	//
	pt->natom =  0;
	pt->njobs = -1;
	pt->matom =  0;
	//
	ok   = -1;
	fp   = fopen ( inpA, "r" );
	if (fp){
		while(NULL!=fgets(sp, BUFSIZE, fp)){
			for(i=0; i<BUFSIZE; ++i) sq[i] = '\0';
			for(i=0; i<strlen(sp);++i) sq[i]=(char)tolower(sp[i]);
			//
			if      (strstr(sq, "microiteration")        ){ pt->njobs = 0; ok = 1; break;}
			else if (strstr(sq, "force calculation")     ){ pt->njobs = 1; ok = 1; break;}
			else if (strstr(sq, "frequency calculation") ){ pt->njobs = 2; ok = 1; break;}
		}
		//
		fclose(fp);
	}
		//
		//
		ok =-1;
	bp   = fopen ( inpB, "r" );
	if (bp){
		ok = 0;
		fclose(bp);
	}
		fp   = fopen ( inpA, "r" );
		bp   = fopen ( inpB, "w" );
			ij   =0;
		while(NULL!=fgets(sp, BUFSIZE, fp)){
			for(i=0; i<BUFSIZE; ++i) sq[i] = '\0';
			for(i=0; i<strlen(sp);++i) sq[i]=(char)tolower(sp[i]);
			//
			j  = strlen(sp);	
			if ( j == 1 ) ++ij;
			if (ij == 2 ) {
				fgets   ( sp,BUFSIZE, fp);
				while(NULL!=fgets(sp, BUFSIZE, fp)){
					j  = strlen(sp);	
					//
					if ( j == 1 ){
						break;
					}
					else         {
						j = sscanf (sp, "%s%d%lf%lf%lf%s%s%d", p1, &p2, &p3, &p4, &p5, p6, p7, &p8);
						if      ( j == 8 ){
							fprintf(bp, "%s\t%i\t%12.7f\t%12.7f\t%12.7f\t%s\t %i\n"
							, p1, p2, p3, p4, p5, p6, p8);
						}
						else if ( j == 6 ){
							p8 = -1;
							fprintf(bp, "%s\t%i\t%12.7f\t%12.7f\t%12.7f\t%s\t %i\n"
							, p1, p2, p3, p4, p5, p6, p8);
						}
						++pt->natom;
					}
				}
			}
		}
		fclose(bp);
		fclose(fp);
	//
	free   (inpB);
	free   (inpA);

	return ok;
}
// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
double **double_Matrix(const int size_a, int size_b){

	int i;
	double **matrix;

	matrix = (double **)malloc(sizeof(double *)*size_a);
	if (!matrix){printf("Allocation Error \n"); exit(1);}
	matrix[0] = (double *)malloc(sizeof(double *)*size_a*size_b);
	if (!matrix){printf("Allocation Error \n"); exit(1);}
	for (i=1;i<size_a;++i){
		matrix[i] = matrix[0] + i*size_b;
	}
	return matrix;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
void  free_double_Matrix(double **matrix){
	free(matrix[0]); free(matrix);
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
double *double_Vector(const int size_a) {

	double *v;
	int a=size_a+1;

	v = (double *)malloc((size_t)a*sizeof(double));
        if (!v) exit(1);
	return v;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
void free_double_Vector(double *v) {
	free ( (char *) (v));
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
void *char_Matrix(int size_a, int size_b){

        char **matrix;
        int i;

        matrix = (char **)malloc(sizeof(char *)*size_a);
        if (!matrix){printf("Allocation Error \n"); exit(1);}
        matrix[0] = (char *)malloc(sizeof(char *)*size_a*size_b);
        if (!matrix){printf("Allocation Error \n"); exit(1);}
        for ( i=1;i<size_a;++i){
                matrix[i] = matrix[0] + i*size_b;
        }
        return matrix;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
void  free_char_Matrix(char **matrix){
        free(matrix[0]); free(matrix);
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
int *int_Vector(const int n)
{
	int *v;
	int a=n+1;

	v = (int *)malloc((size_t)a*sizeof(int));
	if (!v) exit(1);

	return v;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
void free_int_Vector(int *v)
{
    free ( (char *) (v));
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
int **int_Matrix(const int size_a, int size_b){

	int i;
	int **matrix;

	matrix = (int **)malloc(sizeof(int *)*size_a);
	if (!matrix){printf("Allocation Error \n"); exit(1);}
	matrix[0] = (int *)malloc(sizeof(int *)*size_a*size_b);
	if (!matrix){printf("Allocation Error \n"); exit(1);}
	for (i=1;i<size_a;++i){
		matrix[i] = matrix[0] + i*size_b;
	}
	return matrix;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
void  free_int_Matrix(int **matrix){
	free(matrix[0]); free(matrix);
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
double ***double_3D_malloc(int size_a, int size_b, int size_c){

        int i,j;
        double ***matrix;

        matrix= (double***) malloc( sizeof(double**) * size_a );
        if (!matrix){printf("Allocation Error \n"); exit(1);}
        for ( i= 0; i < size_a; i++ ) {
                matrix[i]= (double**) malloc( sizeof(double*) * size_b );
        }
        if (!matrix){printf("Allocation Error \n"); exit(1);}
                matrix[0][0]= (double*) malloc( sizeof(double) * size_a * size_b * size_c );
                for ( i= 0; i < size_a; i++ ) {
                        for(j=0;j<size_b;j++){
                                matrix[i][j]= (double*) ( matrix[0][0] + ( j + i*size_b ) * size_c );
                        }
                }
        return matrix;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
void double_3D_dealloc(double ***matrix){

        free(matrix[0][0]);
        free(matrix[0]);
        free(matrix);
}
// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
static int     MsMonioM_ReL   ( MsMlib *pt){

	
	FILE   *fp, *bp;
	char   *out, *in, s[BUFSIZE], ss[BUFSIZE], ja[50], jb[50];
	int     mo, ia, len, i, j, k, *index;;

	ok   = -1;

if      ( strstr(pt->MsMmM, "xtb")   ){
	in   = (char*)malloc((3 + 4 + 1)*sizeof(char));
	sprintf(in, "ReL.inp\0");

	bp   = fopen (in, "r");
	if  (bp){

		len  = strlen(pt->inp);
		out  = (char*)malloc((len + 8 + 1)*sizeof(char));
		sprintf(out, "%s_ReL.inp\0", pt->inp);
        
		fp   = fopen (out, "w");
		if  (fp){
			while(NULL!=fgets(s, BUFSIZE, bp)){
				fprintf (fp, "%s", s );
			}
			fclose(fp);
		}
		free(out);
		//
		fclose(bp);
	}
	free(in);
	//
		len  = strlen(pt->inp);
		out  = (char*)malloc((len + 12 + 1)*sizeof(char));
		sprintf(out, "%s_ReL.xyz\0", pt->inp);
        
		fp   = fopen (out, "w");
		if  (fp){
				fprintf (fp, "%i\n", pt->natom);
				fprintf (fp, "\n"            );
			for  ( i=0; i<pt->natom; ++i ){
				fprintf(fp,"%s\t %12.7f\t %12.7f\t %12.7f\n" 
				, pt->Typ[i], pt->XYZ[i][0], pt->XYZ[i][1], pt->XYZ[i][2]);
			}
			fclose(fp);
		}
		free(out);

}
else if ( strstr(pt->MsMmM, "g16")   ){
// 
	ok   = -1;
	mo   = -1;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.chk\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		mo = 1;
		fclose(fp);
	}
	free   (out);
	//

	in   = (char*)malloc((3 + 4 + 1)*sizeof(char));
	sprintf(in, "ReL.inp\0");

	bp   = fopen (in, "r");
	if  (bp){


	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.com\0", pt->inp);

	fp   = fopen (out, "w");
	if  (fp){
							fprintf(fp,"%cchk=%s_ReL\n", '%', pt->inp);
		while(NULL!=fgets(ss, BUFSIZE, bp)){
			for(i=0; i<BUFSIZE   ;++i) s[i] = '\0';
			for(i=0; i<strlen(ss);++i) s[i]=(char)tolower(ss[i]);
			//
			if      ( strchr(s, '#' )){
				for( i = 0; i < strlen(s); ++i ){

					if ( s[i] == '@' && s[i+1] == 't' && s[i+2] == 'a' && s[i+3] == 's' && s[i+4] == 'k' ){
							i = i + 4;
						if      (  pt->njobs == 1  ) fprintf(fp, "Force Nosymm " );
						else if (  pt->njobs == 2  ) fprintf(fp, "Freq  Nosymm " );
						else                         fprintf(fp, "      Nosymm " );
						if      (  mo  == 1 )        fprintf(fp, "Guess=TCheck " );
					}
					else                                                                                  {
							fprintf(fp, "%c", ss[i] );
					}
				}
			}
			else if ( strstr(s, "@@@" )){
				for  ( i=0; i<pt->natom; ++i ){
					if ( pt->Opt[i] !=-1 ){
						fprintf(fp,"%s\t  0 \t%12.7f\t%12.7f\t%12.7f\n" , 
						pt->Typ[i], pt->XYZ[i][0], pt->XYZ[i][1], pt->XYZ[i][2]);
					}
					else                  {
						fprintf(fp,"%s\t -1 \t%12.7f\t%12.7f\t%12.7f\n" , 
						pt->Typ[i], pt->XYZ[i][0], pt->XYZ[i][1], pt->XYZ[i][2]);
					}
				}
			}
			else                        {
							fprintf(fp, "%s", ss   );
			}
		}
		fclose(fp);
	}
	free(out);

		fclose(bp);
	}
	free(in);

}
else if ( strstr(pt->MsMmM, "dftb+") ){
//   DFTB+

	ok   = -1;

}

//
	return ok;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
static  int     get_MsMonioM_ReL  ( MsMlib *pt){

	
	FILE   *fp;
	char   *out, s[BUFSIZE], *a;
	int     len, i, j, idim, ndim;
	double  gx, gy, gz;

	// parameters for Hessian matrix
	FILE   *bp;
	char   *in;
	double  drv;
	//

if      ( strstr(pt->MsMmM, "xtb")   ){
	ok   = -101;
	//
	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"finished run on")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	ok   = -101;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"total energy")){ 
				a = strchr(s,'y');
				pt->En[1] = 0.0;
				sscanf(++a, "%lf", &pt->En[1]);
				ok =   0; 
				break;
			}
		}
		fclose(fp);
	}
	free(out);

	if   ( ok != 0 ) return ok;

	if   ( pt->njobs == 1 || pt->njobs == 2){
	ok   = -101;


	len  = strlen(pt->inp);
	out  = (char*)malloc((      8 + 1)*sizeof(char));
	sprintf(out, "gradient\0"       );
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"energy")){ 

				for ( i = 0; i < pt->natom; ++i){
					fgets(s, BUFSIZE, fp);
				}
				//
				for ( i = 0; i < pt->natom; ++i){
					fgets (s, BUFSIZE, fp);
					sscanf(s, "%lf%lf%lf", &gx, &gy, &gz); 
					if    (pt->Opt[i] == -1 ){
						gx = 0.0;
						gy = 0.0;
						gz = 0.0;
					}
					pt->Grd1[i][0] =    gx;
					pt->Grd1[i][1] =    gy;
					pt->Grd1[i][2] =    gz;
				}
				ok   =   0; 
				system ( "rm gradient" );
			}
		}
		fclose(fp);
	}
	free(out);
	}

	if   ( pt->njobs == 2 ){
	//   get hess //

	ok   = -101;

	len  = strlen(pt->inp);
	out  = (char*)malloc((      7 + 1)*sizeof(char));
	sprintf(out, "hessian\0"       );
	fp   = fopen (out, "r");
	if  (fp){

				len  = strlen(pt->inp);	
				in   = (char*)malloc((len + 8 + 1)*sizeof(char));
				sprintf(in, "%s_ReL.hes\0", pt->inp);
				bp   = fopen (in, "w");

		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"$hessian")){ 

				for ( i = 0; i < pt->natom*3; ++i){
				for ( j = 0; j < pt->natom*3; ++j){
					fscanf (fp, "%lf", &drv);
					if ( j <= i ) fprintf(bp, "%23.17f\n", drv);
				}
				}
				//
				ok   =   0; 
			}
		}
					fclose(bp);
				free(in);
		fclose(fp);
		system ( "rm hessian" );
	}
	free(out);
	}

}
else if ( strstr(pt->MsMmM, "g16")   ){
// Gaussian 
	ok   = -201;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"Normal termination of Gaussian")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	ok   = -201;
	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 9 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.fchk\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"Total Energy")){ 
					a=strchr(s, 'R');
					pt->En[1] = 0.0;
					sscanf(++a, "%lf", &pt->En[1]);
			}
			else if (strstr(s,"Cartesian Gradient")){ 

				for ( i = 0; i < pt->natom; ++i){
					fscanf (fp, "%lf", &gx); 
					fscanf (fp, "%lf", &gy); 
					fscanf (fp, "%lf", &gz); 
					if ( pt->Opt[i] == -1 ){
						gx = 0.0;
						gy = 0.0;
						gz = 0.0;
					}
					pt->Grd1[i][0] =    gx;
					pt->Grd1[i][1] =    gy;
					pt->Grd1[i][2] =    gz;

				}

				ok   =   0; 
			}
			else if (strstr(s,"Cartesian Force Constants")){ 
				ndim =  (pt->natom*3)*(pt->natom*3+1)/2;

				len  = strlen(pt->inp);	
				in   = (char*)malloc((len + 8 + 1)*sizeof(char));
				sprintf(in, "%s_ReL.hes\0", pt->inp);
				bp   = fopen (in, "w");
				if  (bp){

					for  (idim=0; idim < ndim; ++idim){
						fscanf (fp, "%lf"  , &drv); 
						fprintf(bp, "%14.9f\n",  drv);
					}
					fclose(bp);
				}
				free(in);
			}
		}
		fclose(fp);
	}
	free(out);

}
else if ( strstr(pt->MsMmM, "dftb+") ){
// DFTB

	ok   = -151;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_ReL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"DFTB+ running times")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	//   get engrads //
	ok   = -151;
	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 17 + 1)*sizeof(char));
	sprintf(out, "%s_detailed.out.ReL\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"Total Mermin free energy")){ 
				a = strchr(s,':');
				pt->En[1] = 0.0;
				sscanf(++a, "%lf", &pt->En[1]);
				ok =   0; 
				break;
			}
		}
		fclose(fp);
	}
	free(out);

	if   ( ok != 0 ) return ok;

	if   (  pt->njobs == 1 ||  pt->njobs == 2 ){
	ok   = -151;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 17 + 1)*sizeof(char));
	sprintf(out, "%s_detailed.out.ReL\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"Total Forces")){ 

				//
				for ( i = 0; i < pt->natom; ++i){
					fgets (s, BUFSIZE, fp);
					sscanf(s, "%d%lf%lf%lf", &j, &gx, &gy, &gz); 
					if    (pt->Opt[i] == -1 ){
						gx = 0.0;
						gy = 0.0;
						gz = 0.0;
					}
					pt->Grd1[i][0] =   -gx;
					pt->Grd1[i][1] =   -gy;
					pt->Grd1[i][2] =   -gz;
				}

				ok   =   0; 
			}
		}
		fclose(fp);
	}
	free(out);

	}

	if   (  pt->njobs == 2 ){
	//   get hess //

	ok   = -101;

        len  = strlen(pt->inp);
        out  = (char*)malloc(( len + 8 + 4 + 1)*sizeof(char));
        sprintf(out, "%s_hessian.ReL\0", pt->inp);

	fp   = fopen (out, "r");
	if  (fp){
				len  = strlen(pt->inp);	
				in   = (char*)malloc((len + 8 + 1)*sizeof(char));
				sprintf(in, "%s_ReL.hes\0", pt->inp);
				bp   = fopen (in, "w");

				for ( i = 0; i < pt->natom*3; ++i){
				for ( j = 0; j < pt->natom*3; ++j){
					fscanf (fp, "%lf", &drv);
					if ( j <= i ) fprintf(bp, "%23.17f\n", drv);
				}
				}
				//
				ok   =   0; 
					fclose(bp);
				free(in);

		fclose(fp);
	}
	free(out);
	}
}

	return ok;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static int     exe_MsMonioM_ReL  (MsMlib *pt){
		FILE   *fp;
		char   *out;
		int     len;

	ok = -1;

	len  = strlen(pt->inp);

	out  = (char*)malloc((len + 7 + 1)*sizeof(char));
	sprintf(out, "%s_go.csh\0", pt->inp);
	fp   = fopen (out, "w");
	if  (fp){

		fprintf(fp,"%c!/bin/csh -f\n", '#');
		fprintf(fp,"\n");
		fprintf(fp,"\n");
		fprintf(fp,"if (-f %s ) rm -rf %s\n", pt->inp, pt->inp);
		fprintf(fp,"if (-d %s ) rm -rf %s\n", pt->inp, pt->inp);
		fprintf(fp,"mkdir -p %s\n", pt->inp);

if      ( strstr(pt->MsMmM, "xtb")   ){
// XTB
		//
		fprintf(fp,"cp -f %s_ReL.%c ./%s/\n", pt->inp, '*', pt->inp);
		fprintf(fp,"cd ./%s/\n", pt->inp);
		if       ( pt->njobs == 2  ){
		fprintf(fp,"${execxtb} %s_ReL.xyz  --hess --grad  --input %s_ReL.inp > %s_ReL.log \n", pt->inp, pt->inp, pt->inp);
		fprintf(fp,"cp %s_ReL.log    ../\n", pt->inp);
		fprintf(fp,"cp gradient      ../\n"       );
		fprintf(fp,"cp hessian       ../\n"       );
		}
		else if  ( pt->njobs == 1  ){
		fprintf(fp,"${execxtb} %s_ReL.xyz  --grad --input %s_ReL.inp > %s_ReL.log \n", pt->inp, pt->inp, pt->inp);
		fprintf(fp,"cp %s_ReL.log    ../\n", pt->inp);
		fprintf(fp,"cp gradient      ../\n"       );
		}

}
else if ( strstr(pt->MsMmM, "g16")   ){
// Gaussian 

		fprintf(fp,"setenv execgau g16\n");
		fprintf(fp,"setenv execchk formchk\n");
		fprintf(fp,"\n");

		//
		fprintf(fp,"cp -f %s_ReL.%c ./%s/\n", pt->inp, '*', pt->inp);
		fprintf(fp,"cd ./%s/\n", pt->inp);
		fprintf(fp,"${execgau} %s_ReL.com > %s_ReL.log\n", pt->inp, pt->inp);
		fprintf(fp,"${execchk} %s_ReL.chk > /dev/null\n", pt->inp);
		fprintf(fp,"cp %s_ReL.chk        ../\n", pt->inp);
		fprintf(fp,"cp %s_ReL.com        ../\n", pt->inp);
		fprintf(fp,"cp %s_ReL.log        ../\n", pt->inp);
		fprintf(fp,"cp %s_ReL.fchk       ../\n", pt->inp);
		fprintf(fp,"cd ../\n");
		//
		//
}
else if ( strstr(pt->MsMmM, "dftb+") ){
// DFTB
		//
		fprintf(fp,"cp -f %s_dftb_in.hsd ./%s/dftb_in.hsd\n", pt->inp, pt->inp);
		fprintf(fp,"cp -f %c.skf ./%s/\n", '*', pt->inp );
		fprintf(fp,"cp -f %s_charges.bin.ReL ./%s/charges.bin\n", pt->inp, pt->inp );
		fprintf(fp,"cp -f %s_ReL.xyz ./%s/\n", pt->inp, pt->inp );
		fprintf(fp,"cd ./%s/\n", pt->inp);
		fprintf(fp,"${execdftb} > %s_ReL.log \n", pt->inp);
		fprintf(fp,"cp detailed.out     ../%s_detailed.out.ReL\n", pt->inp );
		fprintf(fp,"cp -f charges.bin   ../%s_charges.bin.ReL\n", pt->inp );
		fprintf(fp,"cp -f %s_ReL.log   ../\n", pt->inp );
		if ( pt->njobs == 2 ) fprintf(fp,"cp -f hessian.out   ../%s_hessian.ReL\n", pt->inp );
		//
}

		//
		fprintf(fp,"cd ../\n");
		fprintf(fp,"\n");
		fprintf(fp,"rm -rf  %s\n", pt->inp);
		fprintf(fp,"\n");
		fprintf(fp,"\n");
		fclose (fp);
	}

	free(out);
	
	out  = (char*)malloc((9 + len + 7 + 1)*sizeof(char));
	sprintf(out, "chmod +x %s_go.csh\0", pt->inp);
	system (out);
	free   (out);
	
	out  = (char*)malloc((2 + len + 7 + 1)*sizeof(char));
	sprintf(out, "./%s_go.csh\0", pt->inp);
	system (out);
	free   (out);

	return ok;

}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static int      MsMonioM_MoL  ( MsMlib *pt){

	FILE   *fp, *bp;
	char   *out, *in, s[BUFSIZE], ss[BUFSIZE], ja[50], jb[50];
	int     mo, ia, len, i, j, k, *index;;

	ok   = -1;

if      ( strstr(pt->MsMmM, "xtb")   ){
// XTB 
	in   = (char*)malloc((3 + 4 + 1)*sizeof(char));
	sprintf(in, "MoL.inp\0");

	bp   = fopen (in, "r");
	if  (bp){

		len  = strlen(pt->inp);
		out  = (char*)malloc((len + 8 + 1)*sizeof(char));
		sprintf(out, "%s_MoL.inp\0", pt->inp);
        
		fp   = fopen (out, "w");
		if  (fp){
			while(NULL!=fgets(s, BUFSIZE, bp)){
				fprintf (fp, "%s", s );
			}
        
			fclose(fp);
		}
		free(out);
		fclose(bp);
	}
	free(in);


		len  = strlen(pt->inp);
		out  = (char*)malloc((len + 12 + 1)*sizeof(char));
		sprintf(out, "%s_MoL.xyz\0", pt->inp);
        
		fp   = fopen (out, "w");
		if  (fp){
        
				fprintf (fp, "%i\n", pt->matom);
				fprintf (fp, "\n"            );
			for  ( i=0; i<pt->natom; ++i ){
				if      ( strstr(pt->HML[i], "H") || ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){
					fprintf(fp,"%s\t %12.7f\t %12.7f\t %12.7f\n" 
					, pt->Mod[i], pt->MSQ[i][0], pt->MSQ[i][1], pt->MSQ[i][2]);
				}
			}
			fclose(fp);
		}
		free(out);

}
else if ( strstr(pt->MsMmM, "g16")   ){
// Gaussian


	ok   = -1;
	mo   = -1;
	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.chk\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		mo = 1;
		fclose(fp);
	}
	free   (out);
	//

	in   = (char*)malloc(( 3  + 4 + 1)*sizeof(char));
	sprintf(in, "MoL.inp\0");

	bp   = fopen (in, "r");
	if  (bp){


	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.com\0", pt->inp);

	fp   = fopen (out, "w");
	if  (fp){

				fprintf(fp,"%cchk=%s_MoL\n", '%', pt->inp);

		while(NULL!=fgets(ss, BUFSIZE, bp)){
			for(i=0; i<BUFSIZE   ;++i) s[i] = '\0';
			for(i=0; i<strlen(ss);++i) s[i]=(char)tolower(ss[i]);
			//
			if      ( strchr(s, '#' )){
				for( i = 0; i < strlen(s); ++i ){

					if ( s[i] == '@' && s[i+1] == 't' && s[i+2] == 'a' && s[i+3] == 's' && s[i+4] == 'k' ){
							i = i + 4;
					if      (  pt->njobs == 1  ) fprintf(fp, "Force Nosymm " );
					else if (  pt->njobs == 2  ) fprintf(fp, "Freq  Nosymm " );
					else                         fprintf(fp, "      Nosymm " );
					if      (  mo  == 1 )        fprintf(fp, "Guess=TCheck " );
					}
					else                                                                                  {
							fprintf(fp, "%c", ss[i] );
					}
				}
			}
			else if ( strstr(s, "@@@" )){
				for  ( i=0; i<pt->natom; ++i ){
					if      ( strstr(pt->HML[i], "H")                     ){  
						fprintf(fp,"%s\t  0 \t%12.7f\t%12.7f\t%12.7f\n" , 
						pt->Mod[i], pt->MSQ[i][0], pt->MSQ[i][1], pt->MSQ[i][2]);
					}
					else if ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ){  // For Two-Layer
						fprintf(fp,"%s\t  0 \t%12.7f\t%12.7f\t%12.7f\n" , 
						pt->Mod[i], pt->MSQ[i][0], pt->MSQ[i][1], pt->MSQ[i][2]);
					}
				}
			}
			else                        {
							fprintf(fp, "%s", ss   );
			}
		}
		fclose(fp);
	}
	free(out);

		fclose(bp);
	}
	free(in);


}
else if ( strstr(pt->MsMmM, "dftb+") ){
// DFTB+


	ok   = -1;
}
	//

	return ok;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static  int     get_MsMonioM_MoL  ( MsMlib *pt){

	FILE   *fp;
	char   *out, s[BUFSIZE], *a;
	int     len, i, j, k;
	double  gx, gy, gz;

	// parameters for Hessian matrix
	FILE   *bp;
	char   *in;
	double  drv, *dij, zero;
	int     ndim, ii, jj, ij, ji;
	//


if      ( strstr(pt->MsMmM, "xtb")   ){
	ok   = -102;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"finished run on")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	//   get engrads //
	ok   = -102;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"total energy")){ 
				a = strchr(s,'y');
				pt->En[2] = 0.0;
				sscanf(++a, "%lf", &pt->En[2]);
				ok =   0; 
				break;
			}
		}
		fclose(fp);
	}
	free(out);

	if   ( ok != 0 ) return ok;

	if   ( pt->njobs == 1 || pt->njobs == 2 ){
	ok   = -102;

	len  = strlen(pt->inp);
	out  = (char*)malloc((      8 + 1)*sizeof(char));
	sprintf(out, "gradient\0"       );
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"energy")){ 

				for ( i = 0; i < pt->matom; ++i){
					fgets(s, BUFSIZE, fp);
				}
				//
				for ( i = 0; i < pt->natom; ++i){
					if      ( strstr(pt->HML[i], "H")  ||  ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){
					fgets(s, BUFSIZE, fp);
					sscanf (s, "%lf%lf%lf", &gx, &gy, &gz); 
						if      ( pt->Opt[i] == -1 ){
							gx = 0.0;
							gy = 0.0;
							gz = 0.0;
						}
					pt->Grd2[i][0] =    gx;
					pt->Grd2[i][1] =    gy;
					pt->Grd2[i][2] =    gz;
					}
				}

				ok   =   0; 

				system ( "rm gradient" );
			}
		}
		fclose(fp);
	}
	free(out);

	}


	if   ( pt->njobs == 2 ){

	ok   = -102;

	ndim = (pt->matom*3)*(pt->matom*3+1)/2;
	dij  = double_Vector(ndim);

	len  = strlen(pt->inp);
	out  = (char*)malloc((      7 + 1)*sizeof(char));
	sprintf(out, "hessian\0"       );
	fp   = fopen (out, "r");
	if  (fp){

		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"$hessian")){ 
					k =  0;
				for ( i = 0; i < pt->matom*3; ++i){
				for ( j = 0; j < pt->matom*3; ++j){
					fscanf (fp, "%lf", &drv);
					if ( j <= i ) {
						dij[k] = 0.0;
						dij[k] = drv;
						++k;
					}
				}
				}
				//
				ok   =   0; 
			}
		}
		fclose(fp);
	}
	free(out);
	system ( "rm hessian" );

	zero = 0.0;
	len  = strlen(pt->inp);	
	in   = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(in, "%s_MoL.hes\0", pt->inp);
	bp   = fopen (in, "w");
	if   ( bp ){
			k =  0;
		for  (  ii=0; ii < pt->natom; ++ii){
		for  (  ij=0; ij < 3     ; ++ij){

		for  (  jj=0; jj < pt->natom; ++jj){
		for  (  ji=0; ji < 3     ; ++ji){


			if ( jj*3 + ji <= ii*3 + ij ){

				if ( strstr(pt->HML[ii], "H")  ||  ( strstr(pt->HML[ii], "L") && pt->Con[ii] != -1 ) ){  
				if ( strstr(pt->HML[jj], "H")  ||  ( strstr(pt->HML[jj], "L") && pt->Con[jj] != -1 ) ){ 
					fprintf(bp, "%23.17f\n",  dij[k]);
					++k;
				}
				else {
					fprintf(bp, "%23.17f\n",  zero);
				}
				}
				else{
					fprintf(bp, "%23.17f\n",  zero);
				}
			}

		}
		}


		}
		}
			fclose(bp);
			ok  = 0;
	}
	free(in);
	free_double_Vector(dij);
	
	}
}
else if ( strstr(pt->MsMmM, "g16")   ){

// Gaussian
	ok   = -202;
	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"Normal termination of Gaussian")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	//   get engrads //

	ok   = -202;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 9 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.fchk\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"Total Energy")){ 
					a=strchr(s, 'R');
					pt->En[2] = 0.0;
					sscanf(++a, "%lf", &pt->En[2]);
			}
			else if (strstr(s,"Cartesian Gradient")){ 

				for ( i = 0; i < pt->natom; ++i){
			if      ( strstr(pt->HML[i], "H")  ||  ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){  // For Two-Layer
					fscanf (fp, "%lf", &gx); 
					fscanf (fp, "%lf", &gy); 
					fscanf (fp, "%lf", &gz); 

					if ( pt->Opt[i] == -1 ){
						gx = 0.0;
						gy = 0.0;
						gz = 0.0;
					}
					pt->Grd2[i][0] =    gx;
					pt->Grd2[i][1] =    gy;
					pt->Grd2[i][2] =    gz;

			}
				}

				ok =   0; 
			}
			else if (strstr(s,"Cartesian Force Constants")){ 

				len  = strlen(pt->inp);	
				in   = (char*)malloc((len + 8 + 1)*sizeof(char));
				sprintf(in, "%s_MoL.hes\0", pt->inp);
				bp   = fopen (in, "w");
				if  (bp){

					for  (  ii=0; ii < pt->natom; ++ii){
					for  (  ij=0; ij < 3     ; ++ij){

					for  (  jj=0; jj < pt->natom; ++jj){
					for  (  ji=0; ji < 3     ; ++ji){

						if ( jj*3 + ji <= ii*3 + ij ){

							if      ( strstr(pt->HML[ii], "H")  ||  
								( strstr(pt->HML[ii], "L") && pt->Con[ii] != -1 ) ){  // For Two-Layer
							if      ( strstr(pt->HML[jj], "H")  ||  
								( strstr(pt->HML[jj], "L") && pt->Con[jj] != -1 ) ){  // For Two-Layer
								fscanf (fp, "%lf"  , &drv); 
								fprintf(bp, "%14.9f\n",  drv);
							}
							else {
								drv = 0.0;
								fprintf(bp, "%14.9f\n",  drv);
							}
							}
							else{
								drv = 0.0;
								fprintf(bp, "%14.9f\n",  drv);
							}
						}

					}
					}


					}
					}
					fclose(bp);
				}
				free(in);
			}
		}
		fclose(fp);
	}
	free(out);

}
else if ( strstr(pt->MsMmM, "dftb+") ){

// DFTB
	ok   = -152;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoL.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"DFTB+ running times")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	//   get engrads //
	ok   = -152;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 17 + 1)*sizeof(char));
	sprintf(out, "%s_detailed.out.MoL\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"Total Mermin free energy")){ 
				a = strchr(s,':');
				pt->En[2] = 0.0;
				sscanf(++a, "%lf", &pt->En[2]);
				ok =   0; 
				break;
			}
		}
		fclose(fp);
	}
	free(out);

	if   ( ok != 0 ) return ok;

	if   (  pt->njobs == 1 ||  pt->njobs == 2 ){
	ok   = -152;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 17 + 1)*sizeof(char));
	sprintf(out, "%s_detailed.out.MoL\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"Total Forces")){ 
				//
				for ( i = 0; i < pt->natom; ++i){
					if      ( strstr(pt->HML[i], "H")  ||  
						( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){
						fgets(s, BUFSIZE, fp);
						sscanf (s, "%d%lf%lf%lf", &j, &gx, &gy, &gz); 
							if      ( pt->Opt[i] == -1 ){
								gx = 0.0;
								gy = 0.0;
								gz = 0.0;
							}
						pt->Grd2[i][0] =   -gx;
						pt->Grd2[i][1] =   -gy;
						pt->Grd2[i][2] =   -gz;
					}
				}
				ok   =   0; 
			}
		}
		fclose(fp);
	}
	free(out);
	}
	if   ( pt->njobs == 2 ){

	ok   = -102;

	ndim = (pt->matom*3)*(pt->matom*3+1)/2;
	dij  = double_Vector(ndim);

	len  = strlen(pt->inp);
        out  = (char*)malloc(( len + 8 + 4 + 1)*sizeof(char));
        sprintf(out, "%s_hessian.MoL\0", pt->inp);

	fp   = fopen (out, "r");
	if  (fp){

					k =  0;
				for ( i = 0; i < pt->matom*3; ++i){
				for ( j = 0; j < pt->matom*3; ++j){
					fscanf (fp, "%lf", &drv);
					if ( j <= i ) {
						dij[k] = 0.0;
						dij[k] = drv;
						++k;
					}
				}
				}
				ok   =   0; 
		fclose(fp);
	}
	free(out);

	zero = 0.0;
	len  = strlen(pt->inp);	
	in   = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(in, "%s_MoL.hes\0", pt->inp);
	bp   = fopen (in, "w");
	if   ( bp ){
			k =  0;
		for  (  ii=0; ii < pt->natom; ++ii){
		for  (  ij=0; ij < 3     ; ++ij){

		for  (  jj=0; jj < pt->natom; ++jj){
		for  (  ji=0; ji < 3     ; ++ji){


			if ( jj*3 + ji <= ii*3 + ij ){

				if ( strstr(pt->HML[ii], "H")  ||  ( strstr(pt->HML[ii], "L") && pt->Con[ii] != -1 ) ){  
				if ( strstr(pt->HML[jj], "H")  ||  ( strstr(pt->HML[jj], "L") && pt->Con[jj] != -1 ) ){ 
					fprintf(bp, "%23.17f\n",  dij[k]);
					++k;
				}
				else {
					fprintf(bp, "%23.17f\n",  zero);
				}
				}
				else{
					fprintf(bp, "%23.17f\n",  zero);
				}
			}

		}
		}


		}
		}
			fclose(bp);
			ok  = 0;
	}
	free(in);
	free_double_Vector(dij);
	
	}
}

	return ok;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static  int     exe_MsMonioM_MoL  ( MsMlib *pt){

		FILE   *fp;
		char   *out;
		int     len;

	ok   = -1;
	len  = strlen(pt->inp);

	out  = (char*)malloc((len + 7 + 1)*sizeof(char));
	sprintf(out, "%s_go.csh\0", pt->inp);
	fp   = fopen (out, "w");
	if  (fp){

		fprintf(fp,"%c!/bin/csh -f\n", '#');
		fprintf(fp,"\n");
		fprintf(fp,"\n");
		fprintf(fp,"if (-f %s ) rm -rf %s\n", pt->inp, pt->inp);
		fprintf(fp,"if (-d %s ) rm -rf %s\n", pt->inp, pt->inp);
		fprintf(fp,"mkdir -p %s\n", pt->inp);


if      ( strstr(pt->MsMmM, "xtb")   ){
// XTB 
		//
		fprintf(fp,"cp -f %s_MoL.%c ./%s/\n", pt->inp, '*', pt->inp);
		fprintf(fp,"cd ./%s/\n", pt->inp);
		if       ( pt->njobs == 2 ){
		fprintf(fp,"${execxtb} %s_MoL.xyz  --hess --grad  --input %s_MoL.inp > %s_MoL.log \n", pt->inp, pt->inp, pt->inp);
		fprintf(fp,"cp %s_MoL.log    ../\n", pt->inp);
		fprintf(fp,"cp gradient      ../\n"       );
		fprintf(fp,"cp hessian       ../\n"       );
		}
		else if  ( pt->njobs == 1 ){
		fprintf(fp,"${execxtb} %s_MoL.xyz  --grad --input %s_MoL.inp > %s_MoL.log \n", pt->inp, pt->inp, pt->inp);
		fprintf(fp,"cp %s_MoL.log    ../\n", pt->inp);
		fprintf(fp,"cp gradient      ../\n"       );
		}
			//
}
else if ( strstr(pt->MsMmM, "g16")   ){
// Gaussian 

		fprintf(fp,"setenv execgau g16\n");
		fprintf(fp,"setenv execchk formchk\n");
		fprintf(fp,"\n");

		//
		fprintf(fp,"cp -f %s_MoL.%c ./%s/\n", pt->inp, '*', pt->inp);
		fprintf(fp,"cd ./%s/\n", pt->inp);
		fprintf(fp,"${execgau} %s_MoL.com > %s_MoL.log\n", pt->inp, pt->inp);
		fprintf(fp,"${execchk} %s_MoL.chk > /dev/null\n", pt->inp);
		fprintf(fp,"cp %s_MoL.chk        ../\n", pt->inp);
		fprintf(fp,"cp %s_MoL.com        ../\n", pt->inp);
		fprintf(fp,"cp %s_MoL.log        ../\n", pt->inp);
		fprintf(fp,"cp %s_MoL.fchk       ../\n", pt->inp);
		fprintf(fp,"cd ../\n");
		//
		//
}
else if ( strstr(pt->MsMmM, "dftb+") ){
// DFTB
		//
		fprintf(fp,"cp -f %s_dftb_in.hsd ./%s/dftb_in.hsd\n", pt->inp, pt->inp);
		fprintf(fp,"cp -f %c.skf ./%s/\n", '*', pt->inp );
		fprintf(fp,"cp -f %s_charges.bin.MoL ./%s/charges.bin\n", pt->inp, pt->inp );
		fprintf(fp,"cp -f %s_MoL.xyz ./%s/\n", pt->inp, pt->inp );
		fprintf(fp,"cd ./%s/\n", pt->inp);
		fprintf(fp,"${execdftb} > %s_MoL.log \n", pt->inp);
		fprintf(fp,"cp detailed.out     ../%s_detailed.out.MoL\n", pt->inp );
		fprintf(fp,"cp -f charges.bin   ../%s_charges.bin.MoL\n", pt->inp );
		fprintf(fp,"cp -f %s_MoL.log   ../\n", pt->inp );
		if ( pt->njobs == 2 ) fprintf(fp,"cp -f hessian.out   ../%s_hessian.MoL\n", pt->inp );
		//
}


			fprintf(fp,"cd ../\n");
			fprintf(fp,"\n");
			fprintf(fp,"\n");
			fprintf(fp,"\n");
		fclose (fp);
	}

	free(out);
	
	out  = (char*)malloc((9 + len + 7 + 1)*sizeof(char));
	sprintf(out, "chmod +x %s_go.csh\0", pt->inp);
	system (out);
	free   (out);
	
	out  = (char*)malloc((2 + len + 7 + 1)*sizeof(char));
	sprintf(out, "./%s_go.csh\0", pt->inp);
	system (out);
	free   (out);

	ok = 0;
	return ok;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
double          scale_factor (char Sq, char Sm){

		double factor;

	factor = 0.0;

	if     (Sq=='h' && Sm==' ' )return  0.64300;
	else if(Sq=='h' && Sm=='e' )return  0.64300;
	else if(Sq=='l' && Sm=='i' )return  2.45700;
	else if(Sq=='b' && Sm=='e' )return  1.90900;
	else if(Sq=='b' && Sm==' ' )return  1.58700;
	else if(Sq=='c' && Sm==' ' )return  1.43600;
	else if(Sq=='n' && Sm==' ' )return  1.20900;
	else if(Sq=='o' && Sm==' ' )return  1.09600;
	else if(Sq=='f' && Sm==' ' )return  1.02000;
	else if(Sq=='n' && Sm=='e' )return  0.94500;
	else if(Sq=='n' && Sm=='a' )return  2.98600;
	else if(Sq=='m' && Sm=='g' )return  2.64600;
	else if(Sq=='a' && Sm=='l' )return  2.40000;
	else if(Sq=='s' && Sm=='i' )return  2.19200;
	else if(Sq=='p' && Sm==' ' )return  2.06000;
	else if(Sq=='s' && Sm==' ' )return  1.89000;
	else if(Sq=='c' && Sm=='l' )return  1.79500;
	else if(Sq=='a' && Sm=='r' )return  1.70100;
	else                        return -1.00000; //error

	return factor;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
double		get_factor( char *s_qm, char *s_link, char *s_mm){

	int i, j, len;
	double dum0;
	char  s[5];
	char ss[5];
	double r_qm, r_link, r_mm, factor;

	for (i = 0; i < 5; ++i) s[i] = '\0';
		strcpy (s, s_qm);

		len = strlen(s);
		for (j=0; j<len; ++j){
			s[j] = (char)tolower(s[j]);
		}

		if      ( len == 1) dum0=scale_factor(s[0],' ');	
		else                dum0=scale_factor(s[0],s[1]);	

		r_qm = dum0;

	for (i = 0; i < 5; ++i) s[i] = '\0';
	for (i = 0; i < 5; ++i)ss[i] = '\0';


		len = strlen(s_link);
		for(i=0; i<len; ++i){
			if (s_link[i] == '-') break;
			ss[i] = s_link[i];
		}
		
		sscanf(ss, "%s", s);
		len = strlen(s);
		for (j=0; j<len; ++j){
			s[j] = (char)tolower(s[j]);
		}

		if      ( len == 1) dum0=scale_factor(s[0],' ');	
		else if ( len == 2) dum0=scale_factor(s[0],s[1]);	

		r_link = dum0;

	for (i = 0; i < 5; ++i) s[i] = '\0';
		strcpy (s, s_mm);
		len = strlen(s);
		for (j=0; j<len; ++j){
			s[j] = (char)tolower(s[j]);
		}

		if      ( len == 1) dum0=scale_factor(s[0],' ');	
		else if ( len == 2) dum0=scale_factor(s[0],s[1]);	

		r_mm = dum0;

		if ( 0.0 < r_qm && 0.0 < r_link && 0.0 < r_mm){
			factor = ( r_qm + r_link ) / ( r_qm + r_mm );
		}
		else { factor = 0.0;}

	return factor;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static int      prep_model_str ( MsMlib *pt){

		int     ok=-1, i, j;
		double  x, y, z, dr;

		for ( i = 0; i < pt->natom; ++i ){

				pt->MSQ[i][0] = 0.0;
				pt->MSQ[i][1] = 0.0;
				pt->MSQ[i][2] = 0.0;

			if      ( strstr(pt->HML[i], "H")                     ){  
				x  =   pt->XYZ[i][0];
				y  =   pt->XYZ[i][1];
				z  =   pt->XYZ[i][2];

				sscanf (pt->Typ[i], "%s", pt->Mod[i]);

			}
			else if ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ){  // For Two-Layer
				j  =   pt->Con[i] - 1;
				dr =   get_factor(pt->Typ[j], "H", pt->Typ[i]);

				x  =   pt->XYZ[j][0] + dr * ( pt->XYZ[i][0] - pt->XYZ[j][0]);
				y  =   pt->XYZ[j][1] + dr * ( pt->XYZ[i][1] - pt->XYZ[j][1]);
				z  =   pt->XYZ[j][2] + dr * ( pt->XYZ[i][2] - pt->XYZ[j][2]);

				sscanf ("H", "%s", pt->Mod[i]);

			}
				pt->MSQ[i][0] = x;
				pt->MSQ[i][1] = y;
				pt->MSQ[i][2] = z;
			ok = 0;
		}
	return ok;
}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static int      MsMonioM_MoH  ( MsMlib *pt){

	FILE   *fp, *bp;
	char   *out, *in, s[BUFSIZE], ss[BUFSIZE];
	int     len, i, chk;

	ok   = -1;
	chk  = -1;


if      ( strstr(pt->MsMqM, "g16")   ){
// Gaussian
// 
	
	len  = strlen(pt->inp);
	in   = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(in, "%s_MoH.chk\0", pt->inp);

	fp   = fopen (in, "rb");
	if      ( fp ){
		chk = 0;
		fclose ( fp );
	}

	free ( in );


	//
	in   = (char*)malloc((3 + 4 + 1)*sizeof(char));
	sprintf(in, "MoH.inp\0");

	bp   = fopen (in, "r");
	if  (bp){

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.com\0", pt->inp);

	fp   = fopen (out, "w");
	if  (fp){
				fprintf(fp,"%cchk=%s_MoH\n", '%', pt->inp);

		while(NULL!=fgets(ss, BUFSIZE, bp)){
			for(i=0; i<BUFSIZE   ;++i) s[i] = '\0';
			for(i=0; i<strlen(ss);++i) s[i]=(char)tolower(ss[i]);
			//
			if      ( strchr(s, '#' )){
				for( i = 0; i < strlen(s); ++i ){

					if ( s[i] == '@' && s[i+1] == 't' && s[i+2] == 'a' && s[i+3] == 's' && s[i+4] == 'k' ){
							i = i + 4;
					if      ( pt->njobs == 1 ) fprintf(fp, "Force Nosymm " );
					else if ( pt->njobs == 2 ) fprintf(fp, "Freq  Nosymm " );
					else                                  fprintf(fp, "      Nosymm " );
					}
					else                                                                                  {
							fprintf(fp, "%c", ss[i] );
					}
				}
				if ( chk == 0 )         fprintf(fp, "guess=tcheck\n");
							
			}
			else if ( strstr(s, "@@@" )){
				for  ( i=0; i<pt->natom; ++i ){
					if      ( strstr(pt->HML[i], "H")                     ){  
						fprintf(fp,"%s\t  0 \t%12.7f\t%12.7f\t%12.7f\n" , 
						pt->Mod[i], pt->MSQ[i][0], pt->MSQ[i][1], pt->MSQ[i][2]);
					}
					else if ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ){  // For Two-Layer
						fprintf(fp,"%s\t  0 \t%12.7f\t%12.7f\t%12.7f\n" , 
						pt->Mod[i], pt->MSQ[i][0], pt->MSQ[i][1], pt->MSQ[i][2]);
					}
				}
			}
			else                        {
						fprintf(fp, "%s", ss   );
			}
		}
		fclose(fp);
	}
	free(out);

		fclose(bp);
	}
	free(in);

}
else if ( strstr(pt->MsMqM, "orca")   ){
// ORCA 
	//
	in   = (char*)malloc((3 + 4 + 1)*sizeof(char));
	sprintf(in, "MoH.inp\0");

	bp   = fopen (in, "r");
	if  (bp){


	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.inp\0", pt->inp);

	fp   = fopen (out, "w");
	if  (fp){


		while(NULL!=fgets(ss, BUFSIZE, bp)){
			for(i=0; i<BUFSIZE   ;++i) s[i] = '\0';
			for(i=0; i<strlen(ss);++i) s[i]=(char)tolower(ss[i]);
			//
			if      ( strchr(s, '!' )){
				for( i = 0; i < strlen(s); ++i ){

					if ( s[i] == '@' && s[i+1] == 't' && s[i+2] == 'a' && s[i+3] == 's' && s[i+4] == 'k' ){
							i = i + 4;
					if      (  pt->njobs == 1 ) fprintf(fp, "EnGrad   NoUsesym " );
					else if (  pt->njobs == 2 ) fprintf(fp, "Freq EnGrad NoUsesym " );
					else                                  fprintf(fp, "         NoUseSym " );
					}
					else                                                                                  {
							fprintf(fp, "%c", ss[i] );
					}
				}
			}
			else if ( strstr(s, "@@@" )){
				for  ( i=0; i<pt->natom; ++i ){
					if      ( strstr(pt->HML[i], "H") || ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){
						fprintf(fp,"%s\t %12.7f\t %12.7f\t %12.7f\n" , 
						pt->Mod[i], pt->MSQ[i][0], pt->MSQ[i][1], pt->MSQ[i][2]);
					}
				}
			}
			else                        {
							fprintf(fp, "%s", ss   );
			}
		}
		fclose(fp);
	}
	free(out);

		fclose(bp);
	}
	free(in);

}


	return ok;

}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static  int     get_MsMonioM_MoH  ( MsMlib *pt){

	
	FILE   *fp;
	char   *out, s[BUFSIZE], *a;
	int     len, i;
	double  gx, gy, gz;

	// parameters for Hessian matrix
	int     ii, ij, jj, ji, j, k, l, dum, ndim;
	FILE   *bp;
	char   *in;
	double  drv, *dij, zero, v[5];

if      ( strstr(pt->MsMqM, "g16")   ){
	ok   = -203;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"Normal termination of Gaussian")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	printf( "Normal termination of Gaussian\n");

	//   get engrads //

	ok   = -203;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 9 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.fchk\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"Total Energy")){ 
					a=strchr(s, 'R');
					pt->En[3] = 0.0;
					sscanf(++a, "%lf", &pt->En[3]);
			}
			else if (strstr(s,"Cartesian Gradient")){ 

				for ( i = 0; i < pt->natom; ++i){
			if      ( strstr(pt->HML[i], "H")  ||  ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){  // For Two-Layer
					fscanf (fp, "%lf", &gx); 
					fscanf (fp, "%lf", &gy); 
					fscanf (fp, "%lf", &gz); 
					if ( pt->Opt[i] == -1 ){
						gx = 0.0;
						gy = 0.0;
						gz = 0.0;
					}
					pt->Grd3[i][0] =    gx;
					pt->Grd3[i][1] =    gy;
					pt->Grd3[i][2] =    gz;
			}

				}

				ok =   0; 
			}
			else if (strstr(s,"Cartesian Force Constants")){ 

				len  = strlen(pt->inp);	
				in   = (char*)malloc((len + 8 + 1)*sizeof(char));
				sprintf(in, "%s_MoH.hes\0", pt->inp);
				bp   = fopen (in, "w");
				if  (bp){


					for  (  ii=0; ii < pt->natom; ++ii){
					for  (  ij=0; ij < 3     ; ++ij){

					for  (  jj=0; jj < pt->natom; ++jj){
					for  (  ji=0; ji < 3     ; ++ji){

						if ( jj*3 + ji <= ii*3 + ij ){

							if      ( strstr(pt->HML[ii], "H")  ||  
								( strstr(pt->HML[ii], "L") && pt->Con[ii] != -1 ) ){  // For Two-Layer
							if      ( strstr(pt->HML[jj], "H")  ||  
								( strstr(pt->HML[jj], "L") && pt->Con[jj] != -1 ) ){  // For Two-Layer
								fscanf (fp, "%lf"  , &drv); 
								fprintf(bp, "%14.9f\n",  drv);
							}
							else {
								drv = 0.0;
								fprintf(bp, "%14.9f\n",  drv);
							}
							}
							else{
								drv = 0.0;
								fprintf(bp, "%14.9f\n",  drv);
							}
						}

					}
					}
					}
					}
					fclose(bp);
				}
				free(in);
			}
		}
		fclose(fp);
	}
	free(out);
}
else if ( strstr(pt->MsMqM, "orca")   ){
// ORCA

	ok   = -93;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"ORCA TERMINATED NORMALLY")){ ok =   0; break;}
		}
		fclose(fp);
	}
	free(out);

	if      ( ok != 0 ) return ok;

	printf( "Normal termination of Gaussian\n");

	//   get engrads //
	ok   = -93;

	len  = strlen(pt->inp);
	out  = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.log\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if (strstr(s,"FINAL SINGLE POINT ENERGY")){ 
				a = strchr(s,'Y');
				pt->En[3] = 0.0;
				sscanf(++a, "%lf", &pt->En[3]);
				ok =   0; 
				break;
			}
		}
		fclose(fp);
	}
	free(out);

	if   ( ok != 0 ) return ok;

	//   get engrads //

	if   (  pt->njobs == 1 ||  pt->njobs == 11 ||  pt->njobs == 2 ||  pt->njobs == 12){
	ok   = -93;


	len  = strlen(pt->inp);
	out  = (char*)malloc(( len + 11 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.engrad\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){
		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"# The current gradient in Eh/bohr")){ 
					fgets(s, BUFSIZE, fp);
				//
				for ( i = 0; i < pt->natom; ++i){
					if      ( strstr(pt->HML[i], "H")  ||  ( strstr(pt->HML[i], "L") && pt->Con[i] != -1 ) ){
					fgets(s, BUFSIZE, fp);
					sscanf (s, "%lf", &gx); 
					fgets(s, BUFSIZE, fp);
					sscanf (s, "%lf", &gy); 
					fgets(s, BUFSIZE, fp);
					sscanf (s, "%lf", &gz); 
					if      ( pt->Opt[i] == -1 ){
						gx     =    0.0;
						gy     =    0.0;
						gz     =    0.0;
					}
					pt->Grd3[i][0] =    gx;
					pt->Grd3[i][1] =    gy;
					pt->Grd3[i][2] =    gz;
					}
				}

				ok   =   0; 
			}
		}
		fclose( fp);
//		remove(out);
	}
	free(out);

	}

	if   (  pt->njobs == 2 ||  pt->njobs == 12 ){
	//   get hess //

	ok   = -93;

	len  = strlen(pt->inp);
	out  = (char*)malloc(( len + 9 + 1)*sizeof(char));
	sprintf(out, "%s_MoH.hess\0", pt->inp);
	fp   = fopen (out, "r");
	if  (fp){

		ndim =  (pt->matom*3)*(pt->matom*3+1)/2;
		dij  = double_Vector(ndim);

		while(NULL!=fgets(s, BUFSIZE, fp)){
			if      (strstr(s,"$hessian")){ 
						fgets(s, BUFSIZE, fp);
					j = 0;
				while ( j < pt->matom * 3 ){
						fgets(s, BUFSIZE, fp);
					for     ( i = 0; i < pt->matom * 3; ++i ){
						fgets(s, BUFSIZE, fp);
						k = sscanf(s, "%d%lf%lf%lf%lf%lf", &dum, &v[0], &v[1], &v[2], &v[3], &v[4] );
						//
						k--;
						for  ( l = 0; l < k; ++l){
							
							ii = i + 1;
							jj = j + l + 1;
							if ( jj <= ii ){
								ij      =  ii -jj;
								ji      =  ii*(ii+1)/2-ij-1; 
								dij[ji] =  0.0;
								dij[ji] = v[l]; 
							}
						}
					}
					j = j + 5;
				}
				break;
			}
		}
		fclose(fp);
		//
		zero = 0.0;
		len  = strlen(pt->inp);	
		in   = (char*)malloc((len + 8 + 1)*sizeof(char));
		sprintf(in, "%s_MoH.hes\0", pt->inp);
		bp   = fopen (in, "w");
		if   ( bp ){
			k =  0;
			for  (  ii=0; ii < pt->natom; ++ii){
			for  (  ij=0; ij < 3     ; ++ij){

			for  (  jj=0; jj < pt->natom; ++jj){
			for  (  ji=0; ji < 3     ; ++ji){


				if ( jj*3 + ji <= ii*3 + ij ){

					if      ( strstr(pt->HML[ii], "H")  ||  
						( strstr(pt->HML[ii], "L") && pt->Con[ii] != -1 ) ){  // For Two-Layer
					if      ( strstr(pt->HML[jj], "H")  ||  
						( strstr(pt->HML[jj], "L") && pt->Con[jj] != -1 ) ){  // For Two-Layer
						fprintf(bp, "%23.17f\n",  dij[k]);
						++k;
					}
					else {
						fprintf(bp, "%23.17f\n",  zero);
					}
					}
					else{
						fprintf(bp, "%23.17f\n",  zero);
					}
				}

			}
			}

			}
			}
			fclose(bp);
			ok  = 0;
		}
		free_double_Vector(dij);
		free(in);
	}
	remove (out);
	free   (out);
	}


}
	ok = 0;
	return ok;

}
//
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
static int     exe_MsMonioM_MoH  (MsMlib *pt){

		FILE   *fp;
		char   *out;
		int     len;

	ok   = -1;
	len  = strlen(pt->inp);

	out  = (char*)malloc((len + 7 + 1)*sizeof(char));
	sprintf(out, "%s_go.csh\0", pt->inp);
	fp   = fopen (out, "w");
	if  (fp){

		fprintf(fp,"%c!/bin/csh -f\n", '#');
		fprintf(fp,"\n");
		//fprintf(fp,"cd %s/\n", PwD);
		fprintf(fp,"\n");
		fprintf(fp,"if (-f %s ) rm -rf %s\n", pt->inp, pt->inp);
		fprintf(fp,"if (-d %s ) rm -rf %s\n", pt->inp, pt->inp);
		fprintf(fp,"mkdir -p %s\n", pt->inp);
		fprintf(fp,"\n");

if      ( strstr(pt->MsMqM, "g16")   ){
		fprintf(fp,"setenv execgau g16\n");
		fprintf(fp,"setenv execchk formchk\n");
		fprintf(fp,"\n");

// Gaussian
		//
		fprintf(fp,"cp -f %s_MoH.%c ./%s/\n", pt->inp, '*', pt->inp);
		fprintf(fp,"cd ./%s/\n", pt->inp);
		fprintf(fp,"${execgau} %s_MoH.com > %s_MoH.log\n", pt->inp, pt->inp);
		fprintf(fp,"${execchk} %s_MoH.chk > /dev/null\n", pt->inp);
		fprintf(fp,"cp %s_MoH.chk        ../\n", pt->inp);
		fprintf(fp,"cp %s_MoH.com        ../\n", pt->inp);
		fprintf(fp,"cp %s_MoH.log        ../\n", pt->inp);
		fprintf(fp,"cp %s_MoH.fchk       ../\n", pt->inp);
		fprintf(fp,"cd ../\n");
		//
		//
}
else if ( strstr(pt->MsMqM, "orca")   ){
// ORca

		fprintf(fp,"cp -f %s_MoH.%c ./%s/\n", pt->inp, '*', pt->inp);
		fprintf(fp,"cd ./%s/\n", pt->inp);
		fprintf(fp,"${execorca} %s_MoH.inp > %s_MoH.log\n", pt->inp, pt->inp);
		fprintf(fp,"cp %s_MoH.log      ../ >%c /dev/null \n", pt->inp, '&');
		fprintf(fp,"cp %s_MoH.gbw      ../ >%c /dev/null \n", pt->inp, '&');
		fprintf(fp,"cp %s_MoH.ges      ../ >%c /dev/null \n", pt->inp, '&');
		fprintf(fp,"cp %s_MoH.hess     ../ >%c /dev/null \n", pt->inp, '&');
		fprintf(fp,"cp %s_MoH.engrad   ../ >%c /dev/null \n", pt->inp, '&');
		fprintf(fp,"cd ../\n");

}

		fprintf(fp,"\n");
		fprintf(fp,"rm -rf  %s\n", pt->inp);
		fprintf(fp,"\n");
		fprintf(fp,"\n");
		fclose (fp);
	}

	free(out);
	
	out  = (char*)malloc((9 + len + 7 + 1)*sizeof(char));
	sprintf(out, "chmod +x %s_go.csh\0", pt->inp);
	system (out);
	free   (out);
	
	out  = (char*)malloc((2 + len + 7 + 1)*sizeof(char));
	sprintf(out, "./%s_go.csh\0", pt->inp);
	system (out);
	free   (out);

	ok   =  0;
	return ok;
}
// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
static int              engrad_MsMonioM   (MsMlib *pt){

		int     ok=-1;
		int     i, j;
		double  dr;
		double  ax, ay, az;
		double  gx, gy, gz;

				pt->En[0] = pt->En[1]-pt->En[2]+pt->En[3];

		for ( i = 0; i < pt->natom; ++i ){

						gx = 0.0;
						gy = 0.0;
						gz = 0.0;
			if      ( strstr(pt->HML[i],"H") && pt->Opt[i] != -1 ){

						gx = pt->Grd1[i][0] - pt->Grd2[i][0] + pt->Grd3[i][0];
						gy = pt->Grd1[i][1] - pt->Grd2[i][1] + pt->Grd3[i][1];
						gz = pt->Grd1[i][2] - pt->Grd2[i][2] + pt->Grd3[i][2];
						ax = 0.0;
						ay = 0.0;
						az = 0.0;
				for ( j = 0; j < pt->natom; ++j ){

					if ( i == pt->Con[j] - 1 ){
						dr = get_factor(pt->Typ[i], "H" , pt->Typ[j]);
						dr = 1.0 - dr;

						ax = ax - (pt->Grd2[j][0] - pt->Grd3[j][0])*dr;
						ay = ay - (pt->Grd2[j][1] - pt->Grd3[j][1])*dr;
						az = az - (pt->Grd2[j][2] - pt->Grd3[j][2])*dr;
					}
				}
						gx = gx + ax;
						gy = gy + ay;
						gz = gz + az;
			}
			else if ( strstr(pt->HML[i],"H") && pt->Opt[i] != -1 ){
				continue;
			}
			else if ( strstr(pt->HML[i],"L") && pt->Con[i] != -1 ){

						ax = 0.0;
						ay = 0.0;
						az = 0.0;

						j  = pt->Con[i] - 1;

						dr = get_factor(pt->Typ[j], "H", pt->Typ[i]);
						gx = pt->Grd1[i][0] - (pt->Grd2[i][0] - pt->Grd3[i][0])*dr;
						gy = pt->Grd1[i][1] - (pt->Grd2[i][1] - pt->Grd3[i][1])*dr;
						gz = pt->Grd1[i][2] - (pt->Grd2[i][2] - pt->Grd3[i][2])*dr;

						gx = gx + ax;
						gy = gy + ay;
						gz = gz + az;
			}
			else if ( strstr(pt->HML[i],"L") && pt->Con[i] == -1 ){

						gx = pt->Grd1[i][0] - (pt->Grd2[i][0] - pt->Grd3[i][0]);
						gy = pt->Grd1[i][1] - (pt->Grd2[i][1] - pt->Grd3[i][1]);
						gz = pt->Grd1[i][2] - (pt->Grd2[i][2] - pt->Grd3[i][2]);

			}
						pt->Grd0[i][0] = gx;
						pt->Grd0[i][1] = gy;
						pt->Grd0[i][2] = gz;
		}
						ok = 0;
	return ok;
}
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 
static  int                prep_Hessian     ( MsMlib *pt){


	return -1;
}
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 
static int                 hess_OnioM( MsMlib *pt){


	return -1;
}
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 
static  int hessian_max(int a, int b){


	return -1;
}
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 
static  int hessian_min(int a, int b){


	return -1;
}
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 
static  int hessian_index(int a, int b){


	return -1;
}
// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
static int    MsMonioM_out (MsMlib *pt){

	FILE    *fp;
	char    *inpA;
	double   ZZZ;
	int      len;
	int      i, j, ij;

	//
	len  = strlen(pt->inp);
	inpA = (char*)malloc((len + 8 + 1)*sizeof(char));
	sprintf(inpA, "%s_MO.fchk\0", pt->inp);
	//
	ok   = -1;
	fp   = fopen ( inpA, "w" );
	if (fp){
			fprintf ( fp, "Total Energy                               R     %21.17e\n", pt->En[0]);
			fprintf ( fp, "Cartesian Gradient\n" );
		printf ( "Total Energy                               R     %21.17e\n", pt->En[0]);
		printf ( "Cartesian Gradient\n" );
			ij= 0;	
		for  (  i = 0; i < pt->natom; ++i ){
		for  (  j = 0; j < 3        ; ++j ){
			ZZZ = pt->Grd0[i][j];
			fprintf ( fp, " %- 14.8e", ZZZ);
		printf (  " %- 14.8e", ZZZ);
			ok = -1;
			++ij;
			if ( ij == 5 ) {fprintf( fp, "\n" ); ij = 0; printf("\n"); ok = 0;}
		}
		}
		if  ( ok == -1 ){
			fprintf ( fp, "\n" );
			printf (  "\n" );
		}
		//
		printf ( "Cartesian Force Constants\n" );
			fprintf ( fp, "Cartesian Force Constants\n" );
			j = (pt->natom*3)*(pt->natom*3+1)/2;
			ij= 0;	
		for  (  i = 0; i < j; ++i ){
			fprintf ( fp, " %- 14.8e", ZZZ);
		printf (  " %- 14.8e", ZZZ);
			ok = -1;
			++ij;
			if ( ij == 5 ) {fprintf( fp, "\n" ); ij = 0; printf("\n"); ok = 0;}
		}
		if  ( ok == -1 ){
			fprintf ( fp, "\n" );
			printf (  "\n" );
		}
		//
		fclose(fp);
	}
	//
	//
	free   (inpA);
	ok = 0;
	//
	return ok;
}
// 
// ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
// 
