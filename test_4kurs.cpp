#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include  <gsl/gsl_complex_math.h>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <fstream>

using namespace std;

void FillMatrix_M(double** Matr, int k)
{
 int i, j;
 gsl_matrix * m0 = gsl_matrix_alloc (k, k);
 gsl_matrix * m1 = gsl_matrix_alloc (k, k);
 gsl_matrix * m2 = gsl_matrix_alloc (k, k);
 gsl_matrix * m3 = gsl_matrix_alloc (k, k);
 gsl_matrix * m4 = gsl_matrix_alloc (k, k);
 gsl_matrix * m5 = gsl_matrix_alloc (k, k);

 for (i = 0; i < k; i++){
   for (j = 0; j < k; j++){
     if (i==j) {
       gsl_matrix_set (m0, i, j, 1 );
     }
     else {
       gsl_matrix_set (m0, i, j, 0 );
     }
   }
 }
gsl_matrix_scale(m0, 15.0/16.0);
 for (i = 0; i < k; i++){
   for (j = 0; j < k; j++){
     if ((i==j+1) || (i==j-1)) {
       gsl_matrix_set (m1, i, j, 0.5 );
     }
     else {
       gsl_matrix_set (m1, i, j, 0 );
     }
   }
 }
 gsl_matrix_set (m1, 1, 0, 1 );
gsl_matrix_scale(m1, -0.75);
 for (i = 0; i < k; i++){
   for (j = 0; j < k; j++){
     if ((i==j+2) || (i==j-2)) {
       gsl_matrix_set (m2, i, j, 0.5 );
     }
     else {
       gsl_matrix_set (m2, i, j, 0 );
     }
   }
 }
 gsl_matrix_set (m2, 2, 0, 1 );
 gsl_matrix_set (m2, 1, 1, 0.5 );
gsl_matrix_scale(m2, -3.0/16.0);

gsl_matrix_add(m0,m1);
gsl_matrix_add(m0,m2);
 for (i = 0; i < k; i++){
   for (j = 0; j < k; j++){
     Matr[i][j] = gsl_matrix_get(m0,i,j);
   }
 }
}

void FillMatrix_M1(double** Matr, int k){
  gsl_matrix * m0 = gsl_matrix_alloc (k, k);
  for (int i = 0; i < k; i++){
    for (int j = 0; j < k; j++){
      if (i==j) {
        gsl_matrix_set (m0, i, j, -3 );
      }
      else {
        gsl_matrix_set (m0, i, j, 0 );
      }
    }
  }
  for (int i = 0; i < k; i++){
    for (int j = 0; j < k; j++){
      Matr[i][j] = gsl_matrix_get(m0,i,j);
    }
  }
}

void FillMatrix(double** Matr, int k)
{
  for (int i = 0; i < k; i++)
  {
    for (int j = 0; j < k; j++)
    {
      if(j<=i){
        Matr[i][j]=0;
      }
      else{
        if (i==0){
          if (j%2==0){
            Matr[i][j]=0;
          }
          else{
            Matr[i][j]=j;
          }
        }
        if (i==1){
          if (j%2==1){
            Matr[i][j]=0;
          }
          if (j%2==0){
            Matr[i][j]=2*j;
          }

        }
        if (i%2==0&&i>1){
          Matr[i][j]=Matr[0][j]*2;
        }
        if (i%2==1&&i>1){
          Matr[i][j]=Matr[1][j];
        }

      }
    }
  }
}

void Square(double** Matr,int n1){
  double **Matr1=new double* [n1];
  for (int i = 0; i < n1; i++)
    {
        Matr1[i] = new double [n1];
    }
  for (int i = 0; i < n1; i++)
  {
    for (int j = 0; j < n1; j++)
    {
      Matr1[i][j]=0;
     for (int k = 0; k < n1 ; k++){
               Matr1[i][j] += Matr[i][k] * Matr[k][j];
               //cout<<Matr1[i][j]<<endl;
    }
  }
}
for (int i = 0; i < n1; i++)
{
  for (int j = 0; j < n1; j++)
  {
    Matr[i][j]= Matr1[i][j];
  }
}
}

float chebeshev(int n,float x)
{
        //cout <<"im here"<<n<<endl;
       if (n == 0)
       {
             return 1;
       }
       else if (n == 1)
       {
             return x;
       }
       else
       {
             float polinom;
             polinom = 2 * x * chebeshev(n-1, x)-chebeshev(n-2, x);
            // cout<<"polinom "<<polinom<<endl;
             return polinom;
       }
}

int main (void)
{
  //ofstream file("cheb.dat",  ios_base::app);
  ofstream file("kurs2.dat", ios_base::app);
  ofstream file1("kurs3.dat", ios_base::app);
  ofstream file_new("CHEB_NEW.dat");
  ofstream file_new1("CHEB1_NEW.dat");
  int n=16;
  double k_0 = 0.5, Re = 10000, S = 800, ctg =57.3;
  double p;
  //while (k_0<0.1){
  p = 3.0*ctg + k_0*k_0*S*Re;
  float solution = 6.0/5.0*k_0*Re - k_0*p/3.0;
  cout << "= " << 6.0/5.0*k_0*Re - k_0*p/3.0<< endl;
  cout << "P =" << p << endl;
  double U = -3.0, U_0 = 1.5;
    double **D=new double* [n];
    for (int i = 0; i < n; i++)
      {
          D[i] = new double [n];
      }

      double **M=new double* [n];
      for (int i = 0; i < n; i++)
        {
            M[i] = new double [n];
        }
        double **M1=new double* [n];
        for (int i = 0; i < n; i++)
          {
              M1[i] = new double [n];
          }

          double **A=new double* [n];
          for (int i = 0; i < n; i++)
            {
                A[i] = new double [n];
            }
            double **C=new double* [n];
            for (int i = 0; i < n; i++)
              {
                  C[i] = new double [n];
              }

              double **B=new double* [n];
              for (int i = 0; i < n; i++)
                {
                    B[i] = new double [n];
                }

                double **A_new=new double* [n-2];
                for (int i = 0; i < n-2; i++)
                  {
                      A_new[i] = new double [n-2];
                  }

                  double **B_new=new double* [n-2];
                  for (int i = 0; i < n-2; i++)
                    {
                        B_new[i] = new double [n-2];
                    }
                    double **C_new=new double* [n-2];
                    for (int i = 0; i < n-2; i++)
                      {
                          C_new[i] = new double [n-2];
                      }


                      double **A_new2=new double* [2*(n-2)];
                      for (int i = 0; i < 2*(n-2); i++)
                        {
                            A_new2[i] = new double [2*(n-2)];
                        }
                        double **B_new2=new double* [2*(n-2)];
                        for (int i = 0; i < 2*(n-2); i++)
                          {
                              B_new2[i] = new double [2*(n-2)];
                          }
                          double **C_new2=new double* [2*(n-2)];
                          for (int i = 0; i < 2*(n-2); i++)
                            {
                                C_new2[i] = new double [2*(n-2)];
                            }
                            double **L=new double* [2*(n-2)];
                            for (int i = 0; i < 2*(n-2); i++)
                              {
                                  L[i] = new double [2*(n-2)];
                              }
    FillMatrix(D,n);
    FillMatrix_M(M,n);
    FillMatrix_M1(M1,n);
  /*  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        cout << D[i][j] << "\t";
      }
      cout << endl;
    }*/
  /*  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        cout << M1[i][j] << "\t";
      }
      cout << endl;
    }*/
    Square(D, n);
  cout << endl;
    for (int i = 0; i < n-2; i++)
    {
      for (int j = 0; j < n; j++)
      {
        D[i][j]=4.0*D[i][j];
      //  cout << D[i][j] << "\t";
      }
      //cout << endl;
    }

    /*cout <<">>>>>> M <<<<<<<<<" <<endl;
      for (int i = 0; i < n-2; i++)
      {
        for (int j = 0; j < n; j++)
        {
          cout << M[i][j] << "\t";
        }
        cout << endl;
      }*/

      double **MD_2 = new double* [n];
      for (int i = 0; i < n; i++)
      {
        MD_2[i] = new double[n];

        for (int j = 0; j < n; j++)
        {
          //c[i][j] = 0;
          for (int k = 0; k < n; k++)
            MD_2[i][j] += M[i][k] * D[k][j];
        }
      }
      // Вывод матрицы произведения
    /*  cout << ">>>>>>>>>>>>>" << endl;
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          cout << MD_2[i][j] << " ";
        }
        cout << endl;
      }*/

      double **Mk_2=new double* [n];
      for (int i = 0; i < n; i++)
        {
            Mk_2[i] = new double [n];
        }

        for (int i = 0; i < n; i++){
        {
          for (int j = 0; j < n; j++)
             Mk_2[i][j] = M[i][j] * k_0 * k_0;
        }
      }

      for (int i = 0; i < n; i++){

        for (int j = 0; j < n; j++){
           A[i][j] = MD_2[i][j] -  Mk_2[i][j] - M1[i][j];
      }
    }
    for (int i = 0; i < n; i++){

      for (int j = 0; j < n; j++){
         A[i][j] = A[i][j]*k_0*Re;
    }
    }
    /*  cout << ">>>>>>>>>>>>> A <<<<<<<<<<<<<<<<" << endl;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        cout << A[i][j] << " ";
      cout << endl;
    }*/
cout << endl;

  for (int j = 0; j < n; j++){
    A[n-3][j] = k_0*p/U*(4.0*pow(-1,j)*(pow(j,4) - pow(j,2))/3.0 + k_0*k_0*pow(-1,j)) + k_0*2.0*pow(-1,j+1)*j*j*Re*U_0;
    A[n-4][j] = U_0*4.0*pow(-1,j)*(pow(j,4) - pow(j,2))/3.0 + pow(-1,j)*(k_0*k_0*U_0 - U);
  }
  /*cout << ">>>>>>>>>>>>> A after using boundary conditions <<<<<<<<<<<<<<<<" << endl;
  for (int i = 0; i < n; i++)
  {
  for (int j = 0; j < n; j++)
    cout << A[i][j] << " ";
  cout << endl;
}*/

    for (int i = 0; i < n; i++){

      for (int j = 0; j < n; j++){
        C[i][j] = D[i][j];
        if (i == j){
          C[i][j]=D[i][j]-k_0 * k_0;
        }
    }
  }
  for (int i = 0; i < n; i++){

    for (int j = 0; j < n; j++){
       C[i][j] = C[i][j]*k_0*Re;
  }
  }
/*cout << ">>>>>>>>>>>>> C <<<<<<<<<<<<<<<<" << endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      cout << C[i][j] << " ";
    cout << endl;
  }*/
  for (int j = 0; j < n; j++){
    C[n-3][j] = k_0*Re*2.0*pow(-1,j+1)*j*j;
    C[n-4][j] = k_0*k_0*pow(-1,j) + 4.0*pow(-1,j)*(pow(j,4) - pow(j,2))/3.0;
  }
  //cout << ">>>>>>>>>>>>> C after using boundary conditions <<<<<<<<<<<<<<<<" << endl;
  /*  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        cout << C[i][j] << " ";
      cout << endl;
    }*/


cout << endl;
double **kD_2=new double* [n];
for (int i = 0; i < n; i++)
  {
      kD_2[i] = new double [n];
  }

  for (int i = 0; i < n; i++){
  {
    for (int j = 0; j < n; j++)
       kD_2[i][j] = D[i][j] * k_0 * 2.0 * k_0;
  }
  }
  cout << endl;
  FillMatrix(D,n);
  Square(D, n);
  /*for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      cout << D[i][j] << " ";
    cout << endl;
  }*/

  Square(D, n);
for (int i = 0; i < n; i++)
{
  for (int j = 0; j < n; j++){
    D[i][j]=16.0*D[i][j];
    //cout << D[i][j] << " ";
//  cout << endl;
}
}
cout << endl;
for (int i = 0; i < n; i++){

  for (int j = 0; j < n; j++){
     B[i][j] = D[i][j] - kD_2[i][j];
     if (i == j){
       B[i][j]=B[i][j] + k_0 * k_0 * k_0 * k_0;
     }
}
}

cout << ">>>>>>>>>>>>> B <<<<<<<<<<<<<<<<" << endl;
for (int i = 0; i < n; i++)
{
  for (int j = 0; j < n; j++)
    cout << B[i][j] << " ";
  cout << endl;
}
  for (int j = 0; j < n; j++){
     B[n-3][j] = -3.0*k_0*k_0*2.0*pow(-1,j+1)*j*j + pow(-1,j+3)*8.0*j*j*(j*j-1)*(j*j-4)/15.0;
     B[n-4][j] = 0;
}
cout << ">>>>>>>>>>>>> B after using boundary conditions <<<<<<<<<<<<<<<<" << endl;
/*for (int i = 0; i < n; i++)
{
  for (int j = 0; j < n; j++){
    cout << B[i][j] << " ";
  }
  cout << endl;
}*/
cout << ">>>>>>>>>>>>> Using boundary conditions 1,2 <<<<<<<<<<<<<<<<" << endl;
for (int i = 0; i < n-2; i++)
{
  for (int j = 2; j < n; j++){
    A_new[i][j-2] = A[i][0]*(j*j-1) - A[i][1]*(j*j) + A[i][j];
}
}
cout << endl;
/*cout <<" A " <<endl;
cout << endl;
for (int i = 0; i < n-2; i++)
{
  for (int j = 0; j < n-2; j++)
    cout << A_new[i][j] << " ";
  cout << endl;
}*/
cout << endl;
cout << endl;
for (int i = 0; i < n-2; i++)
{
  for (int j = 2; j < n; j++){
    B_new[i][j-2] = B[i][0]*(j*j-1) - B[i][1]*(j*j) + B[i][j];
}
}
cout << endl;
cout << endl;
/*cout <<" B " <<endl;
cout << endl;
cout << endl;
for (int i = 0; i < n-2; i++)
{
  for (int j = 0; j < n-2; j++)
    cout << B_new[i][j] << " ";
  cout << endl;
}*/

for (int i = 0; i < n-2; i++)
{
  for (int j = 2; j < n; j++){
    C_new[i][j-2] =C[i][0]*(j*j-1) - C[i][1]*(j*j) + C[i][j];
}
}
cout << endl;
cout << endl;
cout << endl;
/*cout <<" C " <<endl;
cout << endl;
for (int i = 0; i < n-2; i++)
{
  for (int j = 0; j < n-2; j++)
    cout << C_new[i][j] << " ";
  cout << endl;
}*/

for (int i = 0; i < 2*(n-2); i++)
{
  for (int j = 0; j < 2*(n-2); j++){
    A_new2[i][j] = 0;

  }
}

for (int i = 0; i < 2*(n-2); i++)
{
  for (int j = 0; j < 2*(n-2); j++){
    C_new2[i][j] = 0;

  }
}

for (int i = 0; i < 2*(n-2); i++)
{
  for (int j = 0; j < 2*(n-2); j++){
    B_new2[i][j] = 0;

  }
}

int l = 0;
for (int i = 0; i < n-2; i++)
{
  for (int j = 0; j < n-2; j++){
    A_new2[l][j*2] = A_new[i][j];
    A_new2[l+1][2*j+1] = A_new[i][j];
  }
  l+=2;
  }

cout <<">>>>>" <<endl;
l=0;
/*for (int i = 0; i < 2*(n-2); i++)
{
  for (int j = 0; j < 2*(n-2); j++)
    cout << A_new2[i][j] << " ";
  cout << endl;
}*/

for (int i = 0; i < n-2; i++)
{
  for (int j = 0; j < n-2; j++){
    C_new2[l][j*2] = C_new[i][j];
    C_new2[l+1][2*j+1] = C_new[i][j];
  }
  l+=2;
  }

/*  cout <<">>>>>" <<endl;

  for (int i = 0; i < 2*(n-2); i++)
  {
    for (int j = 0; j < 2*(n-2); j++)
      cout << C_new2[i][j] << " ";
    cout << endl;
  }*/
l=0;
  for (int i = 0; i < n-2; i++)
  {
    for (int j = 0; j < n-2; j++){
      B_new2[l][j*2+1] = -B_new[i][j];
      B_new2[l+1][2*j] = B_new[i][j];
    }
    l+=2;
    }

  /*  cout <<">>>>>" <<endl;

    for (int i = 0; i < 2*(n-2); i++)
    {
      for (int j = 0; j < 2*(n-2); j++)
        cout << B_new2[i][j] << " ";
      cout << endl;
    }
*/
    for (int i = 0; i < 2*(n-2); i++)
    {
      for (int j = 0; j < 2*(n-2); j++){
        L[i][j] = A_new2[i][j] + B_new2[i][j];
      }
      }
    /*  cout <<">>>>> L <<<<<<<<<" <<endl;
      for (int i = 0; i < 2*(n-2); i++)
      {
        for (int j = 0; j < 2*(n-2); j++){
          cout << L[i][j] << " ";
        }
        cout << endl;
      }*/

      gsl_matrix *a = gsl_matrix_alloc (2*(n-2), 2*(n-2));
      gsl_matrix *b = gsl_matrix_alloc (2*(n-2), 2*(n-2));
      double **coeff_real=new double* [n];
        double **coeff_img=new double* [n];
      for (int i = 0; i < n; i++)
        {
            coeff_real[i] = new double [n];
            coeff_img[i] = new double [n];
        }
        double *x= new double [n];

      for (int i = 0; i < 2*(n-2); i++){
        for (int j = 0; j < 2*(n-2); j++)
          gsl_matrix_set (a, i, j, L[i][j]);
        }
          for (int i = 0; i < 2*(n-2); i++){
            for (int j = 0; j < 2*(n-2); j++)
              gsl_matrix_set (b, i, j, C_new2[i][j]);
            }

            gsl_vector_complex *alpha = gsl_vector_complex_alloc (2*(n-2));
            gsl_vector *beta = gsl_vector_alloc (2*(n-2));
            gsl_matrix_complex *evec = gsl_matrix_complex_alloc (2*(n-2), 2*(n-2));
            gsl_vector_complex *eigenValues = gsl_vector_complex_alloc (2*(n-2));
            gsl_eigen_genv_workspace * w =
              gsl_eigen_genv_alloc (2*(n-2));

            gsl_eigen_genv (a,b,alpha,beta, evec, w);

            gsl_eigen_genv_free (w);

            gsl_eigen_genv_sort (alpha,beta, evec,
                                     GSL_EIGEN_SORT_ABS_ASC);
                                     int m =0;
                                     bool flag = true;
                                     for (int i = 0; i < 2*(n-2); i++)
                                       {
                                         flag = true;
                                         double eigenValueReal = GSL_REAL(gsl_complex_div_real(
                                         gsl_vector_complex_get(alpha, i), gsl_vector_get(beta, i)));
                                         gsl_complex eigenValueImg = gsl_complex_div_real(
                                         gsl_vector_complex_get(alpha, i), gsl_vector_get(beta, i));

                                         gsl_vector_complex_view evec_i
                                            = gsl_matrix_complex_column (evec, i);

                                          //  printf ("eigenvalue = %g + %gi\n", GSL_REAL(eigenValueImg), GSL_IMAG(eigenValueImg));

                                            //file<<n<<" eigenvalue = "<<GSL_REAL(eigenValueImg)<<endl;
                                            //  printf ("eigenvector = \n");
                                              int l =0;
                                                for (int j = 0; j <2*(n-2); j+=2)
                                                  {
                                                    gsl_complex z =
                                                      gsl_vector_complex_get(&evec_i.vector, j);
                                                      //cout << GSL_REAL(z)<<endl;

                                                  //  printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
                                                    gsl_complex z1 =
                                                      gsl_vector_complex_get(&evec_i.vector, j+1);
                                                    //printf("%g + %gi\n", GSL_REAL(z1), GSL_IMAG(z1));
                                                    if (GSL_REAL(z)*GSL_IMAG(z1) > 1e-10 && GSL_REAL(z1)*GSL_IMAG(z) < -1e-10) {
                                                      //cout << "-" <<endl;
                                                      GSL_REAL(eigenValueImg) = 0;
                                                      GSL_IMAG(eigenValueImg) = 0;
                                                      flag = false;
                                                  }
                                                  else {
                                                    coeff_real[m][l+2]=(GSL_REAL(z) + (-1)*GSL_IMAG(z1));
                                                    coeff_img[m][l+2]=(GSL_REAL(z1) + GSL_IMAG(z));


                                                    //flag = true;
                                                  //  cout << "m "<<m<<endl;
                                                    //cout <<"l "<<l<<endl;
                                                    l+=1;
                                                  }
                                                  //  printf ("eigenvalue = %g + %gi\n", GSL_REAL(eigenValueImg), GSL_IMAG(eigenValueImg));
                                                  //  coeff_real[i][j+2]=GSL_REAL(z);
                                                  }
                                                    printf ("eigenvalue = %g + %gi\n", GSL_REAL(eigenValueImg), GSL_IMAG(eigenValueImg));
                                                   if (flag == true){
                                                    x[m] = GSL_IMAG(eigenValueImg);
                                                     m+=1;
                                                   }
                                              }
                                              cout<<endl;
                                             for (int i=0; i<n-2; i++){
                                                 for (int j=2; j<n;j++){
                                                   coeff_real[i][0]+=(j*j-1)*coeff_real[i][j];
                                                   coeff_img[i][0]+=(j*j-1)*coeff_img[i][j];
                                                   //j+=2;
                                                   //cout<<i<<" "<<j<<" "<<coeff[i][j]<<endl;
                                                 }
                                               }
                                               for (int i=0; i<n-2; i++){
                                                  for (int j=2; j<n;j++){
                                                    coeff_real[i][1]-=(j*j)*coeff_real[i][j];
                                                    coeff_img[i][1]-=(j*j)*coeff_img[i][j];
                                                  }
                                                }
                                               for (int i=0 ;i<n-2; i++){
                                                  for (int j=0; j<n;j++){
                                                    //coeff[i][0]+=b[j]*coeff[i][j+2];
                                                   // coeff[i][1]+=b[n+j]*coeff[i][j+2];
                                                    //file<<"eigenvector"<<i<<" "<<coeff[i][j]<<endl;
                                                    cout<<i<<" "<<j<<" "<<coeff_real[i][j]<< " "<< coeff_img[i][j]<<endl;
                                                  }
                                                  cout<<endl;
                                                }

                                 cout<<endl;
                                 gsl_vector_complex_free(alpha);
                                 gsl_vector_free(beta);
                                 gsl_matrix_complex_free(evec);
                                 gsl_matrix_free(a);
                                 gsl_matrix_free(b);
                                 for (int i=0 ;i<n-2; i++){

                                      //coeff[i][0]+=b[j]*coeff[i][j+2];
                                     // coeff[i][1]+=b[n+j]*coeff[i][j+2];
                                      //file<<"c_i"<<k_0<<" "<<x[0]<<endl;
                                      cout<<" "<<x[i]<<endl;

                                    cout<<endl;
                                  }
                                  double *alpha_r= new double [n];
                                  double *alpha_i= new double [n];

                                  double *f_r = new double [n];
                                  double *f_i = new double [n];
                                  for (int i=0;i<n-2;i++){
                                    for (int j=0;j<n;j++){
                                  //    cout<<"K"<<k<<endl;
                                      f_r[i] +=coeff_real[i][j] *chebeshev(j,-1);
                                      f_i[i]+=coeff_img[i][j] *chebeshev(j,-1);
                                    //  cout<<"c="<<coeff_img[i][j]<<endl;
                                    }
                                      //cout <<"Cheb " <<cheb[k][i] << endl;
                                    //  cout<<y<<" "<<cheb[k][i]<<endl;
                                  }
                                  for (int i=0;i<n-2;i++){
                                  alpha_r[i] = f_r[i]/(f_r[i]*f_r[i] + f_i[i]*f_i[i]);
                                  alpha_i[i] = -f_i[i]/(f_r[i]*f_r[i] + f_i[i]*f_i[i]);
                                }
                                  cout << "Function "<<f_r[0]<<" "<<f_i[0]<<endl;
                                  cout << "Function "<<f_r[1]<<" "<<f_i[1]<<endl;
                                //  file<<k_0<<" "<<x[0]<<endl;
                                  //file1<<k_0<<" "<<x[1]<<endl;
                                  double t,h;
                                  for (int i=0;i<n-2;i++){
                                    for (int j=0;j<n;j++){
                                  //    cout<<"K"<<k<<endl;
                                  t = coeff_real[i][j];
                                  h = coeff_img[i][j];
                                  coeff_real[i][j]=t*alpha_r[i] - h*alpha_i[i];
                                  coeff_img[i][j]=h*alpha_r[i] + t*alpha_i[i];


                                    //  cout<<"c="<<coeff_img[i][j]<<endl;
                                    }
                                      //cout <<"Cheb " <<cheb[k][i] << endl;
                                    //  cout<<y<<" "<<cheb[k][i]<<endl;
                                  }

                                  float y,phi_1,phi_2,phi_3,phi_0;
                                  int k=0;
                                  int i,j;
                                  double **cheb=new double* [101];
                                  for (int i = 0; i < 101; i++)
                                    {
                                        cheb[i] = new double [n];
                                    }

                                    k=0;
                                    y=0;
                                    while(y<=1){
                                   //cout << chebeshev(5, 2*x-1) << endl;
                                  for (i=0;i<n-2;i++){
                                    for (j=0;j<n;j++){
                                  //    cout<<"K"<<k<<endl;
                                      cheb[k][i] +=coeff_real[i][j]* chebeshev(j, 2*y-1);
                                    //  cout <<"!!!! " <<chebeshev(j, 2*y-1) << endl;
                                    //  cout<<"c="<<coeff_img[i][j]<<endl;
                                    }
                                      //cout <<"Cheb " <<cheb[k][i] << endl;
                                    //  cout<<y<<" "<<cheb[k][i]<<endl;
                                  }
                                  y+=0.01;
                                  k+=1;
                                  //cout<<"x "<<y<<endl;
                                }
                                y=0;
                                float max = -1, max1 = -1;
                        //  for (i=0;i<n-2;i++){
                        for (j=0;j<101;j++){
                          if (abs(cheb[j][1]) > max) {
                            max = abs(cheb[j][1]);
                          }
                          if (abs(cheb[j][0]) > max1) {
                            max1 = abs(cheb[j][0]);
                            cout << max1<<endl;
                          }
                        }
                        cout<<"HELLO MAX "<<max<<endl;
                            for (j=0;j<101;j++){
                              file_new1<<y<<' '<<cheb[j][2]<<endl;
                              file_new<<y<<' '<<cheb[j][4]<<endl;
                              y+=0.01;
                            }
                        //    y=0;
                          //}
                               /* float delta_c, delta_b, delta_d,phi;
                                delta_c = -2.0*solution/3.0 + 1;
                                delta_d = -1.0/9.0*k_0*p + 1.0/2.0*k_0*Re;
                                delta_b = k_0*Re/20.0 - delta_d + 2.0*solution/3.0;
                                phi = -k_0*Re/20.0;
                                cout << phi << ' ' << delta_b <<' '<<delta_c<<' '<<delta_d<<endl;

                              //  k_0+=0.005;
//}*/
                                 return 0;
  }
