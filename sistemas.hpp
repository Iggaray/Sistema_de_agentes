#ifndef SISTEMAS_HPP_INCLUDED
#define SISTEMAS_HPP_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <time.h>

using namespace std;


class Sistema_Multiplicativo
{
    protected:
        int n; //tamaño del sistema
        vector <double> x;
        vector <double> q, lamda, u; //recursos, tasa reseteo, crecimiento, valor reseteo
        vector <double> p_res, mu; //q, lamda escaleadas por tiempo, eficiencia cálculo
        double t, dt; //tiempo de evolucion, paso del sistema
        long seed; //semilla

    public:
        Sistema_Multiplicativo(double q, double lamda, double u, int n);
        Sistema_Multiplicativo(vector <double> q, vector <double> lamda, vector <double> u);

        virtual void paso();

        void transitorio(int pasos);

        int num_agentes();

        vector <double> get_x();

        double xtot();

        double mean();

        double stdev();

        double xmax();

        double entropia();

        double get_t();

        vector<double> top_N(unsigned n);

        int argcum_x(double x);
};


class Replicadores : public Sistema_Multiplicativo
{
    private :
        vector <double> du; //ancho de reseteo [u+- du]

    public:
        Replicadores(double q, double lamda, double u, int n)
        : Sistema_Multiplicativo(q, lamda, u, n)
        {du.resize(n, 0);}

        Replicadores(vector <double> q, vector <double> lamda, vector <double> u, vector <double> du_)
        : Sistema_Multiplicativo(q, lamda, u)
        {du = du_;}

        Replicadores(vector <double> q, vector <double> lamda, vector <double> u)
        : Sistema_Multiplicativo(q, lamda, u)
        {du.resize(q.size(), 0);}

        void paso(); //creo que es lo unico que cambia

        void paso_reset_fijo();

        void paso_sin_reset();
};


#endif // SISTEMAS_HPP_INCLUDED
