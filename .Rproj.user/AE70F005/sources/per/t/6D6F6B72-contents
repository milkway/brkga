#include "RcppArmadillo.h"

const double pi = arma::datum::pi;

//[[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
SEXP Ackleys(const std::vector< double >& X){
  const int N = X.size();
  const arma::vec P = arma::conv_to< arma::vec >::from(X);
  const double subsum1 = arma::sum(P%P);
  const double subsum2 = arma::sum(arma::cos(2*pi*P));
  const double f = -20.*exp(-.2*sqrt(subsum1/(double)N))-exp(subsum2/(double)N)+20.+exp(1); 
  /* note that GAMS uses factor of -0.02 not 0.2  */
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP AluffiPentini(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Aluffi Pentini needs two variables");
  const double f = .25*pow(X[0],4) -.5*X[0]*X[0] + .1*X[0] +.5*X[1]*X[1];
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP BeckerLago(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Becker Lago needs two variables");
  const double f =  pow(fabs(X[0])-5.,2) + pow(fabs(X[1])-5.,2);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Bohachevsky1(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Bohachevsky 1 needs two variables");
  const double f = pow(X[0],2) + 2.*X[1]*X[1] - .3*cos(3*pi*X[0]) -.4*cos(4*pi*X[1]) + .7;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Bohachevsky2(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Bohachevsky 2 needs two variables");
  const double f =  pow(X[0],2) + 2.*X[1]*X[1] - .3*cos(3*pi*X[0])*cos(4*pi*X[1]) + .3;
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP Branin(const std::vector< double >& X,
            const double a=1., 
            const double b=5.1/(4*pi*pi), 
            const double c=5/pi,
            const double d=6., 
            const double e = 10., 
            const double f = 1/(8*pi)
            )
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Branin needs two variables");
  const double r = a*pow((X[1]-b*X[0]*X[0]+c*X[0]-d),2) + e*(1-f)*cos(X[0]) + e;
  return(Rcpp::wrap(r));
}

// [[Rcpp::export]]
SEXP Camel3(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Camel 3 needs two variables");
  const double f = (2 - 1.05*X[0]*X[0] + pow(X[0],4)/6)*X[0]*X[0] + X[0]*X[1] + X[1]*X[1];
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Camel6(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Camel 6 needs two variables");
  const double f = (4.- 2.1*X[0]*X[0] + pow(X[0],4.)/3.)*X[0]*X[0] + X[0]*X[1] + (-4. + 4.*X[1]*X[1])*X[1]*X[1];
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP CosMixN(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  double sum1=0., sum2=0;

  for (auto it = X.begin(); it != X.end(); it++)
  {
    sum1 += cos(5*pi*(*it));
    sum2 += pow((*it),2);
  }
  ;
  const double f =-0.1*sum1 + sum2;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP DekkersAarts(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Dekkers Aarts needs two variables");
  const double f = 100000.*X[0]*X[0] + X[1]*X[1] - pow((X[0]*X[0]+X[1]*X[1]),2) + pow((X[0]*X[0]+X[1]*X[1]),4)/100000.;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Easom(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Easom needs two variables");
  const double f = -cos(X[0])*cos(X[1])*exp(-pow(X[0]-pi,2)-pow(X[1]-pi,2));
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP EMichalewicz(const std::vector< double >& X,
                  const double m = 10.)
  /*Calculate objective function value of x[].*/
{
  const double cost = cos(pi/6.);
  const double sint = sin(pi/6.);
  const int nn = X.size();
  double y[nn];
  int j = 0;
  for (j = 0; j < nn - 1; j += 2)		/* Corrects errors in original
    ICEO test bed */
  {
    y[j] = X[j]*cost - X[j+1]*sint;
    y[j+1] = X[j]*sint + X[j+1]*cost;
  }

  if (j == nn-1)  y[j] = X[j];
  
  double f = 0;
  for (j = 0; j < nn; j++)
    f -= sin(y[j]) * pow(sin((j+1)*y[j]*y[j]/pi),2.0*m);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Expo(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  double s = 0;
  for (auto it = X.begin(); it != X.end(); it++)
  {
    s += (*it)*(*it);
  }
  const double f = -exp(-.5*s);
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP GoldPrice(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("GoldPrice needs two variables");
  double s = 0;
  s = (1+(X[0]+X[1]+1)*(X[0]+X[1]+1)*(19-14*X[0]+3*X[0]*X[0]-14*X[1]+6*X[0]*X[1]+3*X[1]*X[1]));
  const double f = s*(30+(2*X[0]-3*X[1])*(2*X[0]-3*X[1])*(18-32*X[0]+12*X[0]*X[0]+48*X[1]-36*X[0]*X[1]+27*X[1]*X[1]));
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP Griewank(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  double s = 0.;
  double p = 1.;
  for (auto it = X.begin(); it != X.end(); ++it)
  { 
    int index = std::distance(X.begin(), it );
    s += X[index]*X[index];
    p *= cos(X[index]/sqrt((double)(index + 1)));
  }
  const double f = s/4000. - p + 1.;
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP Gulf(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 3) Rcpp::stop("Gulf needs three variables");
  double s = 0., f = 0., u;
  for (int j = 0; j < 99; j++)
  {
    u = 25.+ pow(-50*log(.01*j),.66666);
    s = exp(-pow(u-X[1],X[2])/X[0])-0.01*j;
    f += pow(s,2);
  }
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP Hartman3(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 3) Rcpp::stop("Hartman3 needs three variables");
  int i,j;
  double f = 0;
  double subsum, dist = 0.;
  static double a[4][3] = {
    {3, 10, 30},
    {.1, 10, 35},
    {3, 10, 30},
    {.1, 10, 35}};

  static double c[5] = {1, 1.2, 3, 3.2};
  static double p[4][3] = {
    {.3689, .117, .2673},
    {.4699, .4387, .747},
    {.1091, .8732, .5547},
    {.03815, .5743, .8828}};

  for (i=0; i<5; i++)
  {
    subsum = 0.;
    for (j = 0; j < X.size(); j++)
    {
      subsum = subsum - a[i][j]*pow((X[j]-p[i][j]),2);
    }
    f -= c[i]*exp(subsum);
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Hartman6(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 6) Rcpp::stop("Hartman6 needs six variables");
  int i,j;
  double subsum, dist = 0.;
  static double a[4][6] = {
    {10, 3, 17, 3.5, 1.7, 8},
    {.05, 10, 17, .1, 8, 14},
    {3, 3.5, 1.7, 10, 17, 8},
    {17, 8, .05, 10, .1, 14}};

  static double c[5] = {1, 1.2, 3, 3.2};
  static double p[4][6] = {
    {.1312, .1696, .5569, .0124, .8283, .5886},
    {.2329, .4135, .8307, .3736, .1004, .9991},
    {.2348, .1451, .3522, .2883, .3047, .6650},
    {.4047, .8828, .8732, .5743, .1091, .0381}};

  double f = 0;
  for (i=0; i<5; i++)
  {
    subsum = 0.;
    for (j=0; j < X.size(); j++)
    {
      subsum = subsum - a[i][j]*pow((X[j]-p[i][j]),2);
    }
    f -= c[i]*exp(subsum);
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Hosaki(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Hosaki needs 2 variables");
  double s  = 0;
  s = (1. - 8.*X[0] + 7.*X[0]*X[0] - 7./3*X[0]*X[0]*X[0] + 1./4*X[0]*X[0]*X[0]*X[0]);
  const double f = s*X[1]*X[1]*exp(-X[1]);
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP Kowalik(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("Kowalik needs four variables");
  int i;
  static double a[11] = {.1957, .1947, .1735, .16, .0844, .0627, .0456, .0342, .0323, .0235, .0246};
  static double b[11] = {.25, .5, 1, 2, 4, 6, 8, 10, 12, 14, 16};
  double f = 0.;
  for (i=0; i<11; i++){
    f += pow((a[i]-X[0]*(1+X[1]*b[i])/(1+X[2]*b[i]+X[3]*b[i]*b[i])),2);
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP LM1(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int j, nn;
  double zw1, zw2;
  nn = X.size();
  zw1 = 10*pow(sin(pi*(1+.25*(X[0]+1))),2);
  zw2 = pow(.25*(X[(int)nn-1]+1),2);

  for (j=0; j < (nn-1); j++)
  {
    zw1 += pow(.25*(X[j]+1),2)*(1+pow(sin(pi*(.25*X[j+1])),2));
  }

  const double f = pi/(double)nn*(zw1+zw2);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP McCormic(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("McCormic needs 2 variables");
  const double f = sin((X[0]+X[1])) + pow(X[0]-X[1],2) - 1.5*X[0] + 2.5*X[1] + 1.;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP MeyerRoth(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{  
  if (X.size() != 3) Rcpp::stop("McCormic needs 3 variables");
  int i;
  double num=0., den=0.;
  static double t[5] = {1., 2., 1., 2., .1};
  static double v[5] = {1., 1., 2., 2., 0.};
  static double y[5] = {.126, .219, .076, .126, .186};
  double f = 0.;
  for (i=0; i<5; i++)
  {
    num = X[0]*X[2]*t[i];
    den = 1.+X[0]*t[i]+X[1]*v[i];
    f+= pow(num/den-y[i],2);
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP MieleCantrell(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("McCormic needs 4 variables");
  double sum1 = pow(exp(X[0])-X[1],4) + 100.*pow(X[1]-X[2],6);
  double sum2 = pow(tan(X[2]-X[3]),4) + pow(X[0],8);
  const double f = sum1+sum2;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Modlangerman(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 10) Rcpp::stop("McCormic needs 10 variables");
  int i, j;
  double dx, dist;
  double f = 0;
  static double a[5][10] = {
    {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
    {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
    {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
    {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
    {8.074, 8.777, 3.467, 1.867, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567}};

  static double c[5] = {0.806,0.517,0.1,0.908,0.965};

  for (i=0; i<5; i++)
  {
    dist=0;
    for (j=0; j<X.size(); j++)
    {
      dx=X[j]-a[i][j];
      dist+=dx*dx;
    }
    f-=c[i]*(exp(-dist/pi)*cos(pi*dist));
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP ModRosenbrock(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("ModRosenbrock needs 2 variables");
  const double f =100.*pow((X[1]-X[0]*X[0]),2) + pow((6.4*pow(X[1]-0.5,2) - X[0] - 0.6),2);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP MultiGauss(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("MultiGauss needs 2 variables");
  int i;
  static double a[5] = {.5, 1.2, 1., 1., 1.2};
  static double b[5] = {0., 1., 0., -.5, 0.};
  static double c[5] = {0., 0., -.5, 0., 1.};
  static double d[5] = {.1, .5, .5, .5, .5};

  double f = 0.;
  for (i=0; i<5; i++)
  {
    f -= a[i]*exp(-(pow(X[0]-b[i],2)+pow((X[1]-c[i]),2))/pow(d[i],2));
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Neumaier2(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("Neumaier2 needs 4 variables");
  int i,k;
  double sum1=0.0;
  static double b[4] = {8.0, 18.0, 44.0, 114.0};
  double f = 0.;
  for (k = 0; k < X.size(); k++)
  {
    sum1 = 0.;
    for (i=0; i < X.size(); i++)
    {
      sum1 += pow(X[i], k + 1);
    }
    f += pow(b[k] - sum1, 2);
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Neumaier3(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int i;
  double obj1 = 0.0, obj2 = 0.0;

  for (i = 0; i < X.size(); i++)
  {
    obj1 += pow(X[i]-1.0,2);
  }
  for (i = 1; i < X.size(); i++)
  {
    obj2 += X[i]*X[i-1];
  }

  const double f = obj1 - obj2;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Paviani(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int j;
  double sum1 = 0.0, prod1 = 1.0, prod = 1.0;

  for (j = 0; j < X.size(); j++)
  {
    prod1 = prod1*X[j];
    sum1 = sum1 + pow((double) log(X[j]-2),2) + pow((double) log(10-X[j]),2);
  }

  prod = pow((double) prod1, 0.2);
  const double f = sum1 - prod;
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP Periodic(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Periodic needs 2 variables");
  const double f = 1. + pow(sin(X[0]),2) + pow(sin(X[1]),2) - 0.1*exp(-X[0]*X[0] - X[1]*X[1]);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP PowellQ(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("PowellQ needs 4 variables");
  double s = (X[0]+10.*X[0])*(X[0]+10.*X[0]) + 5.*(X[2]-X[3])*(X[2]-X[3]);
  const double f = s + pow((X[1]-2.*X[2]),4) + 10.*pow((X[0]-X[3]),4);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP PriceTransistor(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 9) Rcpp::stop("PriceTransisto needs 9 variables");
  int k;
  double sumsqr=0.0, alpha, beta;
  static double g[5][4] = {
    {0.485, 0.752, 0.869, 0.982},
    {0.369, 1.254, 0.703, 1.455},
    {5.2095, 10.0677, 22.9274, 20.2153},
    {23.3037, 101.779, 111.461, 191.267},
    {28.5132, 111.8467, 134.3884, 211.4823}};

  for (k=0; k<4; k++)
  {
    alpha = (1.0-X[0]*X[1])*X[2]*(exp(X[4]*(g[0][k]-0.001*g[2][k]*X[6]-0.001*X[7]*g[4][k]))-1.0) - g[4][k] + g[3][k]*X[1];
    beta = (1.0-X[0]*X[1])*X[3]*(exp(X[5]*(g[0][k]-g[1][k]-0.001*g[2][k]*X[6]+g[3][k]*0.001*X[8]))-1.0)- g[4][k]*X[0] + g[3][k];
    sumsqr += alpha*alpha + beta*beta;
  }
  const double f = pow(X[0]*X[2] - X[1]*X[3],2) + sumsqr;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Rastrigin(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{

  int j;
  double f = 0;

  for (j = 0; j < X.size(); j++)
  {
    f += X[j]*X[j]-10.*cos(2.*pi*X[j])+10.;
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Rosenbrock(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int j, nn;
  double a = 0.,b = 0.;
  double f = 0.;
  nn = X.size();
  for (j = 0; j < (nn - 1); j++)
  {
    a = X[j]*X[j]-X[j+1];
    b=1.-X[j];
    f += 100.*a*a+b*b;
  }
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Salomon(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int j;
  double s = 0;

  for (j = 0; j < X.size(); j++)
  {
    s += X[j]*X[j];
  }
  const double f = -cos(2.*pi*sqrt(s))+.1*sqrt(s)+1.;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Schaffer1(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Schaffer1 needs 2 variables");
  const double num = pow((sin(sqrt(X[0]*X[0]+X[1]*X[1]))),2) - 0.5;
  const double den = pow((1+.001*(X[0]*X[0]+X[1]*X[1])),2);
  const double f  = 0.5 + num/den;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Schaffer2(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 2) Rcpp::stop("Schaffer2 needs 2 variables");
  const double prod1 = pow(X[0]*X[0]+X[1]*X[1], 0.25);
  const double prod2 = pow(50*(X[0]*X[0]+X[1]*X[1]), 0.1);
  const double f = prod1*(sin(sin(prod2)) + 1.0);
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Schubert(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int i,j;
  double prod = 1.;
  for (i = 0; i < X.size(); i++)
  {
    double s = 0.;
    for (j=1; j<=5; j++)
      s += j*cos((j+1)*X[i]+j);
    prod = prod*s;
  }
  return(Rcpp::wrap(prod));
}

// [[Rcpp::export]]
SEXP Schwefel(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  int j;
  double s = 0.;

  for (j=0; j < X.size(); j++)
  {
    s += X[j]*sin(sqrt(fabs(X[j])));
  }
  const double f = -s;
  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Shekel10(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("Shekel10 needs 4 variables");
  int i,j;
  double den;
  static double a[10][4] = {
    {4, 4, 4, 4},{1, 1, 1, 1},
    {8, 8, 8, 8},{6, 6, 6, 6},
    {3, 7, 3, 7},{2, 9, 2, 9},
    {5, 5, 3, 3},{8, 1, 8, 1},
    {6, 2, 6, 2},{7, 3.6, 7, 3.6}};
  static double c[10]={.1,.2,.2,.4,.4,.6,.3,.7,.5,.5};
  double s = 0.;
  for (i=0; i<10; i++) {
    den = 0.;
    for (j = 0; j < X.size(); j++) {
      den = den + pow((X[j]-a[i][j]),2);
    }
    s -=  1.0/(den + c[i]);
  }
  return(Rcpp::wrap(s));
}

// [[Rcpp::export]]
SEXP Shekel5(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("Shekel15 needs 4 variables");
  int i,j;
  double den;
  static double a[5][4] = {{4, 4, 4, 4},{1, 1, 1, 1},{8, 8, 8, 8},{6, 6, 6, 6},{3, 7, 3, 7}};
  static double c[5]={0.1,0.2,0.2,0.4,0.4};
  double s = 0.;
  for (i = 0; i < 5; i++) {
    den = 0.;
    for (j = 0; j< X.size(); j++) {
      den = den + pow((X[j] - a[i][j]), 2);
    }
    s -= 1.0/(den + c[i]);
  }
  return(Rcpp::wrap(s));
}

// [[Rcpp::export]]
SEXP Shekel7(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("Shekel7 needs 4 variables");
  int i,j;
  double den;
  static double a[7][4] = {
    {4, 4, 4, 4},{1, 1, 1, 1},
    {8, 8, 8, 8},{6, 6, 6, 6},
    {3, 7, 3, 7},{2, 9, 2, 9},
    {5, 5, 3, 3}};
  static double c[7]={.1,.2,.2,.4,.4,.6,.3};
  double s = 0.;
  for (i = 0; i < 7; i++) {
    den = 0.;
    for (j = 0; j < X.size(); j++) {
      den = den + pow((X[j]-a[i][j]),2);
    }
    s -= 1.0/(den + c[i]);
  }
  return(Rcpp::wrap(s));
}

// [[Rcpp::export]]
SEXP Shekelfox5(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 10) Rcpp::stop("Shekelfox5 needs 10 variables");
  int i,j;
  double sp, h;
  static double a[30][10] = {
    {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
    {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
    {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
    {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
    {8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
    {7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
    {1.256, 3.605, 8.623, 6.905, 4.584, 8.133, 6.071, 6.888, 4.187, 5.448},
    {8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
    {0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
    {7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
    {0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
    {2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
    {8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
    {2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
    {4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
    {8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
    {8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
    {4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
    {2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
    {6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
    {0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
    {5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
    {3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
    {8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
    {1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
    {0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
    {0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
    {4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
    {9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
    {4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699,
     6.500}};

  static double c[30] = {
    0.806,
    0.517,
    0.1,
    0.908,
    0.965,
    0.669,
    0.524,
    0.902,
    0.531,
    0.876,
    0.462,
    0.491,
    0.463,
    0.714,
    0.352,
    0.869,
    0.813,
    0.811,
    0.828,
    0.964,
    0.789,
    0.360,
    0.369,
    0.992,
    0.332,
    0.817,
    0.632,
    0.883,
    0.608,
    0.326};
  double s = 0.;
  for (j = 0; j < 30; j++)
  {
    sp = 0.0;
    for (i = 0; i < X.size(); i++)
    {
      h = X[i]-a[j][i];
      sp += h*h;
    }
    s -= 1.0/(sp+c[j]);
  }
  return(Rcpp::wrap(s));
}

// [[Rcpp::export]]
SEXP Wood(const std::vector< double >& X)
  /*Calculate objective function value of x[].*/
{
  if (X.size() != 4) Rcpp::stop("Wood needs 4 variables");
  const double f = 100*pow((X[1]-X[0]*X[0]),2) + (1-X[0])*(1-X[0])
  +90*pow((X[3]-X[2]*X[2]),2) + (1-X[2])*(1-X[2])
  +10.1*(pow((X[1]-1),2)+pow((X[3]-1),2))
  +19.8*(X[1]-1)*(X[3]-1);

  return(Rcpp::wrap(f));
}

// [[Rcpp::export]]
SEXP Zeldasine(const std::vector< double >& X,
               const double A = 2.5, 
               const double B = 5.
                   )
  /*Calculate objective function value of x[].*/
{
  int j;
  double prod1=1., prod2=1.;
  const double z = pi/6;

  for (j = 0; j < X.size(); j++)
  {
    prod1 *= sin(X[j]-z);
    prod2 *= sin(B*(X[j]-z));
  }
  const double f = -(A*prod1 + prod2);
  return(Rcpp::wrap(f));
}


// [[Rcpp::export]]
SEXP reg_beta(const std::vector< double >& X, const arma::mat& data){
  //Definindo vetor eta; eta = mX*vb (produto de matriz por vetor)
  int n_rows = data.n_rows;
  int n_cols = data.n_cols; // usar o N
  arma::vec x = arma::conv_to<arma::vec>::from(X);
  arma::vec resp_var = data.col(0);
  arma::mat model_mat = arma::join_rows(arma::ones<arma::vec>(n_rows,1), data.cols(1, n_cols - 1));

  arma::vec eta = model_mat*x.rows(0,x.size()-2);
  arma::vec mu = arma::exp(eta)/(arma::exp(eta) + 1);
  double phi = x[x.size()-1];
  arma::vec temp = arma::lgamma(mu*phi) - arma::lgamma((1-mu)*phi) + (mu*phi - 1)%arma::log(resp_var) + (((1-mu)*phi ) - 1)%arma::log(1 - resp_var);

  const double f  = n_rows*lgamma(phi) + arma::sum(temp);
  return(Rcpp::wrap(f));
}

