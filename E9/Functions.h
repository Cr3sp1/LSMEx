#ifndef __Functions__
#define __Functions__

#include <cmath>
#include <vector>
#include <armadillo>

using namespace std;


class Function
{

public:
  virtual double Eval(double) const = 0;
  virtual ~Function(){};
};


class ScalarFunction
{

public:
  virtual double Eval(vector<double> &x) const = 0;
  virtual double Eval(arma::vec &v) const = 0;
  virtual double Eval(arma::mat &m) const = 0;
  virtual ~ScalarFunction(){};
};


class Sin : public Function
{

public:
  Sin(){};
  double Eval(double x) const override { return sin(x); };
  ~Sin(){};
};


class Es2I : public Function
{

public:
  Es2I(){};

  double Eval(double x) const override
  {
    return M_PI * 0.5 * cos(M_PI * x / 2);
  };

  ~Es2I(){};
};

class Es2P : public Function
{

public:
  Es2P(){};

  double Eval(double x) const override
  {
    return 2*(1-x);
  };

  ~Es2P(){};
};


class Es2Q : public Function
{

public:
  Es2Q(){};

  double Eval(double x) const override
  {
    return 1-sqrt(1-x);
  };

  ~Es2Q(){};
};



#endif // __Function__