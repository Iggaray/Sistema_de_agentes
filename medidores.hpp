#ifndef MEDIDORES_HPP_INCLUDED
#define MEDIDORES_HPP_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <forward_list>
#include <array>

using namespace std;

class Estadistico
{
    protected:
        int paso;

    public:
        Estadistico() : paso(0) //constructor
        { }

        virtual void nuevo_dato(double) = 0;

        virtual void reset();

        virtual double get() = 0;

};

class Media : public Estadistico
{
    private:
        double media;

    public:
        Media() : Estadistico() //constructor
        {media = 0;}

        void nuevo_dato(double x);

        void reset();

        double get();
};

class Desvio_Estandar : public Media
{
    private:
        double mu2;

    public:
        Desvio_Estandar() : Media() //constructor
        {mu2 = 0;} //<x> y <x^2> acumulados

        void nuevo_dato(double x);

        double get();

        double get_mean();

        double get_std();

        void reset();
};

class Mu3 : public Desvio_Estandar
{
    private:
        double mu3;

    public:
        Mu3() : Desvio_Estandar()
        {mu3 = 0;}

        void nuevo_dato(double x);
        void nuevo_dato(vector <double> x);

        double get();

        double get_mean();

        double get_std();

        void reset();

        void print();
};

class Histograma
    {
    ///Representa un histograma con espaciamiento lineal
    protected:
        vector <double> bins;
        vector <long> cuentas;
        int paso, nbins;
        double xmin, xmax;

    public:
        Histograma(double x_min, double x_max, int n);

        virtual int __index__(double x) = 0;

        void nuevo_dato(double x);
        void nuevo_dato(vector <double> x);

        void reset();

        vector <double> get_normalizado();

        vector <double> get_densidad();

        vector <double> get_bins();

        vector <long> get_cuentas();

        void print();

        void plot_densidad(FILE* grafico = NULL, bool log = false);

        void save_densidad(char* nombre_archivo);

    };

class HistoLineal: public Histograma
    {
    public:
        HistoLineal(double x_min, double x_max, int n); //constructor

        int __index__(double x);
    };

class HistoLog: public Histograma
{
    protected:
        double razon;

    public:
        HistoLog(double x_min, double x_max, int n); //constructor

        int __index__(double x);
};

class Plano_de_fases
{
    private:
        double xmin, xmax, nx, ymin, ymax, ny;
        vector <double> grilla;
        HistoLog xbins, ybins;

    public:
        Plano_de_fases(double xmin, double xmax, int nx, double ymin, double ymax, int ny);

        vector <int> __index__(double x, double y);

        vector <double> get_xbins()
        {return xbins.get_bins();}

        vector <double> get_ybins()
        {return ybins.get_bins();}

        void nuevo_dato(double valor, double x, double y);

        void nuevo_dato(double valor, int idx, int idy);

        void print();

        void print_xbins();

        void print_ybins();

        void plot();

        void save(char* nombre_archivo);
};

class Trayectoria
{
    private:
        forward_list < pair<double, double> > trayectoria;

    public:
        Trayectoria();

        void nuevo_dato(double t, double x);

        void save(char* nombre_archivo);

        void plot();

        vector <double> get_x();

        vector <double> get_t();
};
#endif // MEDIDORES_HPP_INCLUDED
