#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

void initialiserCholesky(vector<vector<float>> &matrice,vector<vector<float>>&transposee,vector<float> &second,int &dim);
void afficherMatrice(vector<vector<float>>& matrice,int& dim);
void cholesky(vector<vector<float>>& A,int &dim);
void cholesky(vector<vector<float>>& A,int &dim);
vector<vector<float>> transpose(vector<vector<float>>&A,int &dim );
void resolution(vector<vector<float>> &A,vector<vector<float>> &transpose,vector<float>&B,int &dim);

int main(){
        vector<vector<float>> A,transposee;
        vector<float> B;
        int dim = 0;
        cout<<"Résolvons le système d'équations ci-dessous à l'aide de la méthode de Cholesky:"<<endl;
        initialiserCholesky(A,transposee,B,dim);
        cholesky(A,dim);
        //afficherMatrice(A,dim);
        cout<<" La matrice résultante de la méthode de Cholesky:"<<endl;
        afficherMatrice(A,dim);
        cout<<"Transposée de la matrice triangularisée:"<<endl;
        transposee = transpose(A,dim);
        afficherMatrice(transposee,dim);
        resolution(A,transposee,B,dim);
        return 0;
}

void resolution(vector<vector<float>> &A,vector<vector<float>> &transpose,vector<float>&second,int &dim){
//vecteur qui contient les inconnues du système AY = B 
//On initialise 2 vecteurs pour où Y = transposée(A)*X (X vecteur inconnu)
    vector<float> y,x ;
//initialisation des inconnues
    for(int l=0; l<dim;l++){
        x.push_back(0);
        y.push_back(0);
    }
    float somme,somme1;
    for(int i = 0; i<dim; i++ ){
        somme = 0;
        for(int j =0; j<i; j++){
//produit de la matrice A.X
            somme+=(y[j]*A[i][j]);   
        }
//On obtient le résultat grâce à b[i]-somme/elt de la diagonale de A
        y[i] = (second[i] - somme )/A[i][i];
    }

    for(int i = dim-1; i>=0; i-- ){
        somme1 = 0;
        for(int j =i+1; j<dim; j++){

            somme1+=(x[j]*transpose[i][j]);
        }

        x[i] = (y[i] - somme1 )/transpose[i][i];
    }
    cout<<"Le système admet pour solutions :"<<endl;
//affichage des solutions du système d'équations AX =B
    for(int i=0;i<x.size();i++){
        cout<<'\t'<<x[i]<<endl;
    }
}

vector<vector<float>> transpose(vector<vector<float>>&A,int &dim ){
    vector<vector<float>> result;
    for (int i=0 ;i<dim;i++){
        vector<float> colonne;
        for (int j=0;j<dim;j++){
            colonne.push_back(0);
        }
        result.push_back(colonne);
    }
    for (int i=0 ;i<dim;i++){
        for (int j=0;j<dim;j++){
            result[i][j] = A[j][i];
        }
    }
    return result;
}


void cholesky(vector<vector<float>>& A,int &dim){

    for(int i=0;i<dim;i++){
        
        for(int j=0;j<dim; j++){
           if(j<i){
                int k =0;
                float somme =0;
                while(k<j){
                    somme += A[i][k]*A[j][k];
                    k++;
                }
            //cout<<"somme:"<<j<<somme<<endl;
                A[i][j] = 1/A[j][j]*(A[i][k]- somme);
           }else if(j==i){
            int l =0;
            float sommeDiag = 0;
            while(l<i){   
                sommeDiag+= A[i][l]*A[i][l];
                l++;
            }
            
            A[i][i] = A[i][i] - sommeDiag;
            A[i][i] = sqrt(A[i][i]);
           }else{
                A[i][j] = 0; 
           }
        }
    }
}


void afficherMatrice(vector<vector<float>>& matrice,int& dim){
    for (int i=0;i<dim;i++){
       for (int j=0; j<dim;j++){
            cout<<matrice[i][j]<<"\t";
        }
       cout<<endl;
    }
}

void initialiserCholesky(vector<vector<float>> &matrice,vector<vector<float>> &transposee,vector<float> &second,int &dim){
//initialisation d'une variable afin d'ouvrir le fichier
    ifstream fichier("data.txt");
    if(fichier){
//la première ligne du fichier est déstinéé à être la dimension de la matrice
        fichier>>dim;
//Tant que la variable i sera inférieure à la dimension insérée ci-dessus
//nous considérerons les lignes en dessous de la première comme étant les
//coéfficients de la matrice
        for(int i=0;i<dim;i++){
            vector<float> ligne;
            for(int j=0;j<dim;j++){
                float temp(0);
                fichier>>temp;
                ligne.push_back(temp);
            }
            matrice.push_back(ligne);
        }
//le reste du fichier sera la matrice de B du système d'équation AX = B où A est 
//la matrice à triangulariser
        for(int i=0;i<dim;i++){
            float temp = 0;
            fichier>>temp;
            second.push_back(temp);
        }
        for(int i=0;i<dim;i++){
            vector<float> colonne;
            for(int j=0;j<dim;j++){
                colonne.push_back(0);
            }
            transposee.push_back(colonne);
        }
// on affiche le système d'équations  AX = B
        for (int i= 0;i<dim ; i++){
            cout<<"|\t";
            for(int j =0;j<dim; j++){
                cout<<matrice[i][j]<<"x["<<j<<"] +";
                if(j == dim -1)
                    cout<<matrice[i][j]<<"x["<<j<<"] =";
            }
            cout<< second[i]<<endl;
        }
    }
//En cas d'erreur d'ouverture du fichier, on affiche un message d'erreur
else{
        cout<<"Erreur lors de chargement du fichier"<<endl;
    }
}