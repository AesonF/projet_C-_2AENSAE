#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <vector>


class Matrix {
private:
    int n;
    int m;
    float **p;
    void allocArray();
public:
    Matrix(int a, int b); //constructeur d'une matrice de zeros
    Matrix(int a, int b, std::vector<float> &A); //constructeur pour les tests d'operator
    Matrix id(int a); //matrice identite de taille a
    void load(int i, int j, float a); //insertion d'une valeur a en (i,j)
    float load(int i, int j); //observation de la valeur en (i,j)
    int lin(); //nombre de lignes
    int col(); //nombre de colonnes
    void show(); //affichage d'une matrice
    Matrix operator *(Matrix &A);
    void operator =(Matrix &A);
    void operator *=(Matrix&A); //ATTENTION! Il s'agit bien de this*A, pas A*this
    Matrix copy();
    ~Matrix();
};

Matrix quickExp(Matrix &A, int a);

Matrix transpose(Matrix &A);

void LaTeXShow(Matrix &A); //display for copy-pasting to LaTeX matrix environments

//Les trois prochaines fonctions servent uniquement à inverser une matrice
void dilate(Matrix &A, int a, int b);
void transvect(Matrix &A, int i, int j, float lambda);
void rearrange(Matrix &A, int i);
Matrix inverse(Matrix &A);

#endif // MATRIX_H
