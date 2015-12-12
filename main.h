#if !defined(MAIN_H)
/* ========================================================================
   $File: $
   $Date: $
   $Revision: $
   $Creator: Casey Muratori $
   $Notice: (C) Copyright 2015 by Molly Rocket, Inc. All Rights Reserved. $
   ======================================================================== */

#define MAIN_H
#endif

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
using namespace std;

void interp2 ( double* sig, int N );
void decim2( double* sig, int N );
double* conv( double* sig, int N, double* h, int K );
double* concat( const double* tab1, int N, const double* tab2, int M);
void analyse_haar( double** x, int p);
double* synthese_haar( const double* coeff, int p);
void print_tab( const double* tab, int N );
void add_tab( const double* tab1, int N, const double* tab2, int M );
void analyse_97( double* x, int p );
double* synthese_97( double* x, int p );
void analyse_97_lifting( double* x, int p );
double * synthese_97_lifting( double* x, int p );
void write_signal_to_file( const char* filename, double* x);
double * get_signal_from_file( const char* filename );
double error(const double* tab1, int N, const double* tab2, int M);
double * copy( const double * tab, int N );
int getIdx( int id, int N );
void amr(double* x,int p,int niveau);
void iamr(double* x,int p,int niveau);
void analyse2D_97( double * m, int p, int size);
void synthese2D_97( double* m, int p, int size);
void amr2D_97(double* m,int p,int j);
void iamr2D_97( double* m, int p, int j);
int ecrit_bmp256(const char* fichier, uint32_t largeur, uint32_t hauteur, double* m);
double* charge_bmp256(const char* fichier, uint32_t largeur, uint32_t hauteur);
double mean_band( double* m , int p);
double var_band( double* m, int p, double mean);
double* compute_means( double* m, int p, int j );
double* compute_vars( double* m, int p, int j, double* means);
double* extract_band( double* m, int p, int size);
