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
    Matrix copy(); //"deep copy" d'une matrice (copie les valeurs)
    float load(int i, int j); //observation de la valeur en (i,j)
    int lin(); //nombre de lignes
    int col(); //nombre de colonnes

    void load(int i, int j, float a); //insertion d'une valeur a en (i,j)
    void show(); //affichage d'une matrice

    Matrix operator *(Matrix &A);
    void operator =(Matrix &A);
    void operator *=(Matrix&A); //il s'agit de this*A, pas A*this

    ~Matrix();
};



Matrix quickExp(Matrix &A, int a); //exponentiation rapide
Matrix transpose(Matrix &A); //transposition
void LaTeXShow(Matrix &A); //output texte pour copier-coller dans LaTeX



//Les trois prochaines fonctions servent uniquement à inverser une matrice
void dilate(Matrix &A, int a, int b);
void transvect(Matrix &A, int i, int j, float lambda);
void rearrange(Matrix &A, int i);
Matrix inverse(Matrix &A);

#endif // MATRIX_H
