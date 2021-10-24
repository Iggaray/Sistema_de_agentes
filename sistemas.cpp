#include <iostream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <numeric>
#include <math.h>

#include "sistemas.hpp"
#include "random_generator.hpp"

using namespace std;

///------------ Clase Sistema_Multiplicativo
Sistema_Multiplicativo::Sistema_Multiplicativo(double q_, double lamda_, double u_, int n_)
{///Inicializa un sistema homogeneo
    t = 0;
    dt = 1e-3;
    n = n_;
    seed = time(NULL);
    //inicializo recursos en valor de reseteo
    q.resize(n_, q_);
    lamda.resize(n_, lamda_);
    u.resize(n_, u_);
    x.resize(n_, u_);
    p_res.resize(n_, q_ * dt);
    mu.resize(n_, lamda_ * dt);
}

Sistema_Multiplicativo::Sistema_Multiplicativo(vector <double> q_, vector <double> lamda_, vector <double> u_)
{///Inicializa un sistema general
    t = 0; dt = 1e-3;
    seed = time(NULL);
    q = q_; lamda = lamda_; u = u_;

    x.resize(u.size());
    copy(u.begin(), u.end(), x.begin());

    p_res.resize(q.size());

    mu.resize(lamda.size());

    transform(q.begin(), q.end(), p_res.begin(),
              bind1st(multiplies<double>(), dt));

    transform(lamda.begin(), lamda.end(), mu.begin(),
              bind1st(multiplies<double>(), dt));
}

void Sistema_Multiplicativo::paso()
    {
    for(int i=0; i<n; i++){
        x[i] *= mu[i];
        if (ran2(&seed) < p_res[i]) x[i] = u[i];
        }

    t += dt;
    }

vector <double> Sistema_Multiplicativo::get_x()
    {
    return x;
    }

int Sistema_Multiplicativo::num_agentes()
    {return n;}

void Sistema_Multiplicativo::transitorio(int pasos)
{
    for(int i=0; i<pasos; i++)
        paso();
}

double Sistema_Multiplicativo::xtot()
{
    double suma = 0;
    for(unsigned i=0; i<x.size(); i++)
        suma += x[i];
    return suma;
}

double Sistema_Multiplicativo::mean()
{
    return Sistema_Multiplicativo::xtot() / x.size();
}

double Sistema_Multiplicativo::stdev()
{
    double suma = 0;
    for(unsigned i=0; i<x.size(); i++)
        suma += x[i]*x[i];
    suma -= pow(mean(), 2);
    return sqrt(suma) / x.size();
}

double Sistema_Multiplicativo::xmax()
{
    auto pmax = max_element(x.begin(), x.end());
    return *pmax;
}

vector<double> Sistema_Multiplicativo::top_N(unsigned n)
{
    sort(x.begin(), x.end(), greater<double>());
    if (x.size() < n)
        n = x.size();
    vector<double> top(n);
    for(unsigned i=0; i<n; i++)
        top[i] = x[i];
    return top;
}

int Sistema_Multiplicativo::argcum_x(double frac)
{
    if(frac > 1) {
        printf("Warning: Fracción mayor a uno\n");
        return x.size();
    }

    sort(x.begin(), x.end(), greater<double>()); //ordeno de mayor a menor
    double sum = 0, xtot = frac * accumulate(x.begin(), x.end(), 0.0);
    int n = 0; //contador de agentes
    while(sum < xtot){
        sum += x[n];
        n++;
    }
    return n;
}

double Sistema_Multiplicativo::entropia()
{
    double xT = 0, H = 0;
    for(unsigned i=0; i<x.size(); i++)
        xT += x[i];
    for(unsigned i=0; i<x.size(); i++)
        H -= x[i]/xT * log(x[i]/xT);
    return H;
}

double Sistema_Multiplicativo::get_t()
{
    return t;
}
///------------ Clase Replicadores (general, con reseteos)

void Replicadores::paso()
{
    //fitness promedio * dt
    double bar_f = inner_product(mu.begin(), mu.end(), x.begin(), 0.0);
    //paso determinista + reseteo
    for(unsigned i = 0; i<x.size(); i++){
        x[i] *= 1.0 - bar_f + mu[i];
        if(ran2(&seed) < p_res[i])
            x[i] = u[i] + 2 * du[i] * (ran2(&seed) - 0.5);
    }
    t += dt;
}

void Replicadores::paso_reset_fijo()
{
    //fitness promedio * dt
    double bar_f = inner_product(mu.begin(), mu.end(), x.begin(), 0.0);
    //paso determinista + reseteo
    for(unsigned i = 0; i<x.size(); i++){
        x[i] *= 1.0 - bar_f + mu[i];
        if(ran2(&seed) < p_res[i])
            x[i] = u[i];
    }
    t += dt;
}

void Replicadores::paso_sin_reset()
{
    //sorteo los mu
    for(unsigned i = 0; i < mu.size(); i++){
        mu[i] = lamda[i] * (ran2(&seed)) * dt; //el ancho me lo absorbe la escala temporal
    }
    //fitness promedio * dt
    double bar_f = inner_product(mu.begin(), mu.end(), x.begin(), 0.0);

    //paso determinista
    for(unsigned i = 0; i<x.size(); i++){
        x[i] *= 1.0 - bar_f + mu[i];
    }
    t += dt;
}
