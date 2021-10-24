#include <iostream>
#include <math.h>
#include <algorithm>
#include "medidores.hpp"

using namespace std;

///-------------Clase Estadistico
void Estadistico::reset()
{
    paso = 0;
}

//-------------SubClase Media
void Media::nuevo_dato(double x)
{
    media = (media * paso + x) / (paso + 1);
    paso++;
}

double Media::get()
{
    return media;
}

void Media::reset()
{
    Estadistico::reset();
    media = 0;
}
//-------------SubClase Desvio_Estandar

void Desvio_Estandar::nuevo_dato(double x)
{
    mu2 = (mu2 * paso + x*x) / (paso + 1);
    Media::nuevo_dato(x); //media = (media * paso + x) / (paso + 1);
}

double Desvio_Estandar::get()
{
    double mu = Media::get();
    return sqrt(mu2 - mu*mu);
}

double Desvio_Estandar::get_mean()
{
    return Media::get();
}

void Desvio_Estandar::reset()
{
    Media::reset();
    mu2 = 0;
}

//-------------SubClase Mu3
void Mu3::nuevo_dato(double x)
{
    //Media::nuevo_dato(x); //media = (media * paso + x) / (paso + 1);
    mu3 = (mu3 * paso + x*x*x) / (paso + 1);
    Desvio_Estandar::nuevo_dato(x); //mu2 = (mu2 * paso + x*x) / (paso + 1);
}

void Mu3::nuevo_dato(vector <double> x)
{
    for(unsigned i=0; i<x.size(); i++){
        mu3 = (mu3 * paso + x[i]*x[i]*x[i]) / (paso + 1);
        Desvio_Estandar::nuevo_dato(x[i]); //mu2 = (mu2 * paso + x*x) / (paso + 1);
    }
}

double Mu3::get()
{
    double mu = Media::get();
    double stdev = Desvio_Estandar::get();
    return (mu3 - 3*mu*stdev*stdev - mu*mu*mu) / stdev/stdev/stdev;
}

double Mu3::get_mean()
{
    return Media::get();
}

double Mu3::get_std()
{
    return Desvio_Estandar::get();
}

void Mu3::print()
{
    printf("Media\tStd\tSkewness\n");
    printf("%.5f\t%.5f\t%.5f\n", Media::get(), Desvio_Estandar::get(), get());
}

void Mu3::reset()
{
    mu3 = 0;
    Desvio_Estandar::reset();
}





///--------------Clase Histograma

Histograma::Histograma(double x_min, double x_max, int n) : paso(0), nbins(n), xmin(x_min), xmax(x_max)
    {
    double dx = (x_max - x_min) / (nbins - 1); //espaciamiento lineal

    for(int i=0; i<nbins; i++)
        bins.push_back(x_min + i * dx);

    cuentas.resize(nbins, 0);
    }


void Histograma::nuevo_dato(double x)
{
    cuentas[__index__(x)] += 1;
    paso += 1;
}


void Histograma::nuevo_dato(vector <double> x)
{
    for(unsigned i=0; i<x.size(); i++){
        cuentas[__index__(x[i])] += 1;
        paso += 1;
        }
}


void Histograma::reset()
{
    for(unsigned i=0; i<cuentas.size(); i++)
        cuentas[i] = 0;
    paso = 0;
}


vector <double> Histograma::get_normalizado()
{
    vector <double> cuentas_norm(cuentas.size());
    for(unsigned i=0; i<cuentas.size(); i++)
        cuentas_norm[i] = (double)(cuentas[i]) / paso;
    return cuentas_norm;
}

vector <double> Histograma::get_densidad()
{
    vector <double> dx(bins.size() - 1);
    vector <double> f(bins.size() - 1);
    for(unsigned i=0; i<dx.size(); i++){
        dx[i] = bins[i+1] - bins[i];
        f[i] = ((double)(cuentas[i+1]) + (double)(cuentas[i]))/2 / paso /dx[i];
        }
    return f;
}
void Histograma::print()
{
    for(unsigned i = 0; i<cuentas.size(); i++)
            printf("%.3lf\t%ld\n", bins[i], cuentas[i]);
}

vector <double> Histograma::get_bins()
{
    return bins;
}


vector <long> Histograma::get_cuentas()
{
    return cuentas;
}

void Histograma::plot_densidad(FILE* grafico, bool log)
{
    vector <double> densidad = get_densidad(); //cargo densidad
    if (grafico == NULL)
        grafico = popen( "gnuplot -persist", "w" );
    auto [ymin, ymax] = minmax_element(densidad.begin(), densidad.end());
    if (log)
        fprintf(grafico,"set logscale xy\n");
    fprintf(grafico, "set xrange [%f:%f]\n", xmin, xmax);
    fprintf(grafico, "set yrange [%f:%f]\n", *ymin + 1e-6, *ymax); //esto se puede personalizar para que la chance de pasarse sea minima
    fprintf(grafico,"plot '-' w lp\n");
    for(int i=0; i<nbins-1; i++)
       fprintf(grafico,"%f %f\n", (bins[i] + bins[i+1]) / 2, densidad[i]);
    fprintf(grafico,"e\n");
}

void Histograma::save_densidad(char* nombre_archivo)
{
    vector <double> densidad = get_densidad(); //cargo densidad
    FILE *file = fopen(nombre_archivo, "w" );
    fprintf(file, "bins\tdensity\n");
    for(unsigned i=0; i<densidad.size(); i++)
        fprintf(
            file,
            "%lg\t%lg\n",
            (bins[i+1] + bins[i]) / 2,
            densidad[i]);

    fclose(file);
}


HistoLineal::HistoLineal(double x_min, double x_max, int n)
    : Histograma(x_min, x_max, n) //constructor
        { }

int HistoLineal::__index__(double x)
{
    int index;

    if (x < xmin)
        index = 0;
    else if (x > xmax)
        index = nbins - 1;
    else
        index = (int)((x - xmin) / (xmax - xmin) * (nbins - 1));

    return index;
}

HistoLog::HistoLog(double x_min, double x_max, int n) : Histograma(x_min, x_max, n) //constructor
{
    razon = pow((long double)(xmax/xmin), (long double)(1.0/(nbins-1)));
    double b = 1.0;
    for(int i=0; i<nbins; i++){
        bins[i] = xmin * b;
        b *= razon;
    }
}

int HistoLog::__index__(double x)
{
    int index;

    if (x < xmin)
        index = 0;
    else if (x > xmax)
        index = nbins - 1;
    else
        index = (int)(log(x / xmin) / log(razon));

    return index;
}


///Class Plano_de_fases
Plano_de_fases::Plano_de_fases(double xmin, double xmax, int nx,
                               double ymin, double ymax, int ny)
:xmin(xmin), xmax(xmax), nx(nx), ymin(ymin), ymax(ymax), ny(ny),
 xbins(xmin, xmax, nx), ybins(ymin, ymax, ny)
{grilla.resize(nx * ny);}

void Plano_de_fases::nuevo_dato(double valor, double x, double y)
{
    vector <int> idx_idy = __index__(x, y);
    int idx = idx_idy[0];
    int idy = idx_idy[1];

    grilla[idx + idy * nx] = valor;
}

void Plano_de_fases::nuevo_dato(double valor, int idx, int idy)
{
    grilla[idx + idy * nx] = valor;
}

vector<int> Plano_de_fases::__index__(double x, double y)
{
    int idx = xbins.__index__(x);
    int idy = ybins.__index__(y);

    vector <int> idx_idy = {idx, idy};
    return idx_idy;
}

void Plano_de_fases::print()
{
    printf("\n");
    for(unsigned idx = 0; idx < xbins.get_bins().size(); idx++)
        printf("\t%5.3g", xbins.get_bins()[idx]);
    printf("\n\t");
    for(unsigned i=0; i < xbins.get_bins().size(); i++)
        printf("--------");
    printf("\n");
    for(unsigned idy = 0; idy < ybins.get_bins().size(); idy++){
        printf("%5.3g |\t", ybins.get_bins()[idy]);
        for(unsigned idx = 0; idx < xbins.get_bins().size(); idx++)
            printf("%5.2g\t", grilla[idx + idy * nx]);
        printf("\n%7s\n", "|");
    }
}

void Plano_de_fases::print_xbins()
{
    xbins.print();
}

void Plano_de_fases::print_ybins()
{
    ybins.print();
}

void Plano_de_fases::save(char* nombre_archivo)
{
    FILE* f = fopen(nombre_archivo, "w");
    for(unsigned idy = 0; idy < ybins.get_bins().size(); idy++){
        for(unsigned idx = 0; idx < xbins.get_bins().size(); idx++)
            fprintf(f, "%g,", grilla[idx + idy * nx]);
        fprintf(f, "\n");
    }
    fclose(f);
}

///-------------Clase Trayectoria

Trayectoria::Trayectoria()
{ }

void Trayectoria::nuevo_dato(double t, double x)
{
    trayectoria.emplace_front(t, x);
}

void Trayectoria::plot()
{
    FILE *grafico = popen( "gnuplot -persist", "w" );
    //fprintf(grafico,"set logscale y\n");
    fprintf(grafico,"plot '-' w lp\n");
    for(auto& pos: trayectoria)
        fprintf(grafico,"%lf %lf\n", pos.first, pos.second);
    fprintf(grafico, "set yrange [%f:%f]\n", 1e-7, 1.0);
    fprintf(grafico,"e\n");
}

void Trayectoria::save(char* nom_archivo)
{
    FILE *archivo = fopen(nom_archivo, "w" );
    //fprintf(grafico,"set logscale y\n");
    fprintf(archivo, "t\t x\n");
    for(auto& pos: trayectoria)
        fprintf(archivo, "%lg %lg\n", pos.first, pos.second);
}
