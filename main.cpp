#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <algorithm>

#include "sistemas.hpp"
#include "medidores.hpp"
using namespace std;

void diagrama_de_fases(double q2, int niter)
{
    double lam2 = 1.0;
    vector <double> u = {0.1, 0.1};
    vector <double> du = {0.025, 0.025};

    const double q1_q2_min = 0.1;
    const double q1_q2_max = 100.0;

    const double lam1_lam2_min = 1.0/sqrt(10); //1-0
    const double lam1_lam2_max = sqrt(10.0);

    int nq = 40, nlam = 40;//int nq = 43, nlam = 43;
    Plano_de_fases plano_mu(q1_q2_min, q1_q2_max, nq, lam1_lam2_min, lam1_lam2_max, nlam),
                   plano_mu3(q1_q2_min, q1_q2_max, nq, lam1_lam2_min, lam1_lam2_max, nlam);
    Mu3 mu3_1, mu3_2;

    for(int i=0; i<nq; i++){
        //fijo el q_1 del plano de parametros
        double q1 = plano_mu.get_xbins()[i] * q2;
        vector <double> q = {q1, q2};

        for(int j=0; j<nlam; j++){
            //fijo el lamda_1 del plano de parametros
            double lam1 = plano_mu.get_ybins()[j] * lam2;
            vector <double> lam = {lam1, lam2};

            //Creo sistema de dos replicadores
            Replicadores replicadores(q, lam, u, du);
            replicadores.transitorio(10000);
            mu3_1.reset();
            mu3_2.reset();
            //Simulacion + estadistica
            for(int n=0; n<niter; n++){
                replicadores.paso();

                if(n%500 == 0){
                    mu3_1.nuevo_dato(replicadores.get_x()[0]);
                    mu3_2.nuevo_dato(replicadores.get_x()[1]);
                }
            }
            //Agrego el estadístico al diagrama de fase
            plano_mu.nuevo_dato(mu3_1.get_mean() - mu3_2.get_mean(), i, j);
            plano_mu3.nuevo_dato(mu3_1.get() - mu3_2.get(), i, j);
            printf("\033c"); //limpia pantalla
            plano_mu.print();
        }
    }
    //nomrbo archivo para guardar plano. los bins se supone que los genera el cliente con los
    //datos del nombre
    char* nom_archivo_mu = (char*)malloc(sizeof(char) * 300);
    char* nom_archivo_mu3 = (char*)malloc(sizeof(char) * 300);

    sprintf(
        nom_archivo_mu,
        "/home/nacho/Escritorio/Tesis_Maestria/Simulaciones/Agentes_replicadores/data/replicadores_diagrama_fases_full_mu_q%.1lf_qmin%.2lf_qmax%.2lf_nq%d_lammin%.2lf_lammax%.2lf_nlam%d_u%.2lf_pm%.3lf.csv",
        q2, q1_q2_min, q1_q2_max, nq, lam1_lam2_min, lam1_lam2_max, nlam, u[0], du[0]);

    sprintf(
        nom_archivo_mu3,
        "/home/nacho/Escritorio/Tesis_Maestria/Simulaciones/Agentes_replicadores/data/replicadores_diagrama_fases_full_mu3_q%.1lf_qmin%.2lf_qmax%.2lf_nq%d_lammin%.2lf_lammax%.2lf_nlam%d_u%.2lf_pm%.3lf.csv",
        q2, q1_q2_min, q1_q2_max, nq, lam1_lam2_min, lam1_lam2_max, nlam, u[0], du[0]);

    plano_mu.save(nom_archivo_mu);
    plano_mu3.save(nom_archivo_mu3);

    free(nom_archivo_mu);
    free(nom_archivo_mu3);
}

void simulacion_valor_medio(int n, double* mu, double* sigma, int pasos=(int)(5e6))
{///Simular sistema de n replicadores a lo largo de pasos iteraciones.
    /*Guardar archivo con el promedio temporal del valor medio y de sus
    fluctuaciones temporales*/

    double q = 0.1, lambda = 1.0, u = 0.01;
    //Instancio el sistema de replicadores
    Replicadores sistema(q, lambda, u/n, n);
    sistema.transitorio(10000);

    //Instancio el estadistico que computa mu y sigma
    HistoLog histo_barx(u/n, 1.0/n, 100);
    Desvio_Estandar mu2;

    //Simulo y tomo datos
    for(int i=0; i<pasos; i++){
        sistema.paso_reset_fijo();
        if(i%10 == 0){
            double bar_x = sistema.mean();
            histo_barx.nuevo_dato(bar_x);
            mu2.nuevo_dato(bar_x);
        }
    }

    //nombre histograma
    char* nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../data/replicadores_density_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lambda, u);
    histo_barx.plot_densidad();
    histo_barx.save_densidad(nom_archivo);
    *mu = mu2.get_mean();
    *sigma = mu2.get();
    free(nom_archivo);
}

void simulacion_density_indiv(int n, int pasos=(int)(5e6))
{///Simular sistema de n replicadores a lo largo de pasos iteraciones.
    /*Guardar archivo con el histograma de las densidades individuales*/

    double q = 0.1, lambda = 1.0, u = 0.01;
    //Instancio el sistema de replicadores
    Replicadores sistema(q, lambda, u/n, n);
    sistema.transitorio(10000);

    //Instancio el histograma logarítmico
    HistoLog histo_xi(u/n, 1.0, 100);

    //Simulo y tomo datos
    for(int i=0; i<pasos; i++){
        sistema.paso_reset_fijo();
        if(i%1000 == 0){
            histo_xi.nuevo_dato(sistema.get_x());
            if(i%10000 == 0)
                printf("ETA: ---- %.2lf%%\n", 100.0 * i / pasos);
        }
    }

    //nombre histograma
    char* nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_density_global_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lambda, u);
    histo_xi.plot_densidad();
    histo_xi.save_densidad(nom_archivo);
    free(nom_archivo);
}

void evolucion_repl_global(int n = 1000000, double q=0.1, double lam=1.0, double u=1e-2, int niter=(int)(1e6))
{
    Replicadores sistema(q, lam, u/n, n);
    Trayectoria entropia, xtot, xmax, x2nd, x3rd, n_50prcnt;
    HistoLog distrib_x(u, 1.0, 20);

    //comienzo a simular
    for(int i = 0; i<niter; i++){
        sistema.paso_reset_fijo(); //paso con mismo reseteo para todos
        if (i%2000 == 0){
            xtot.nuevo_dato(sistema.get_t(), sistema.xtot());
            entropia.nuevo_dato(sistema.get_t(), sistema.entropia());
            vector<double>xtop = sistema.top_N(3);
            xmax.nuevo_dato(sistema.get_t(), xtop[0]);
            x2nd.nuevo_dato(sistema.get_t(), xtop[1]);
            x3rd.nuevo_dato(sistema.get_t(), xtop[2]);
            n_50prcnt.nuevo_dato(sistema.get_t(), sistema.argcum_x(0.5));
            if (i%10000 == 0)
                printf("ETA: ------- %.2lf %%\n", 100.0 * i / niter);
            if(fabs(sistema.entropia() - 4.0) < 1e-1)
                distrib_x.nuevo_dato(sistema.get_x());
        }
    }
    distrib_x.plot_densidad(NULL, true);
    distrib_x.save_densidad("../Agentes_replicadores/data/densidad.dat");
    char* nom_archivo = (char*)malloc(sizeof(char) * 300);

    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_xtot_vs_t_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lam, u);
    xtot.save(nom_archivo);
    free(nom_archivo);
    xtot.plot();

    nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_xmax_vs_t_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lam, u);
    xmax.save(nom_archivo);
    free(nom_archivo);
    xmax.plot();

    nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_x2nd_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lam, u);
    x2nd.plot();
    x2nd.save(nom_archivo);
    free(nom_archivo);

    nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_x3rd_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lam, u);
    x3rd.plot();
    x3rd.save(nom_archivo);
    free(nom_archivo);

    nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_entropia_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lam, u);
    entropia.plot();
    entropia.save(nom_archivo);
    free(nom_archivo);

    nom_archivo = (char*)malloc(sizeof(char) * 300);
    sprintf(
        nom_archivo,
        "../Agentes_replicadores/data/replicadores_n50prcnt_N%d_q%.2lf_lam%.2f_u%.2f.csv",
        n, q, lam, u);
    //n_50prcnt.plot();
    n_50prcnt.save(nom_archivo);
    free(nom_archivo);
}



int main(int argc, char** argv)
{
    //simulacion_density_indiv(100, (int)(5e7));
    double q = 0.3, lam = 1.0, u = 1e-2;
    int iter = (int)(1e8);
    for(int i=1; i < 5; i++){
        evolucion_repl_global(3 * pow(10, i), q, lam, u, iter);
    }

/*
    double q2 = 1.0;
    vector <double> u = {0.01, 0.01};
    vector <double> du = {0.0025, 0.0025};
    vector <double> lam = {1.0, 1.0};

    const double q1_q2_min = 1e-1;
    const double q1_q2_max = 1e1;
    const int n_histogramas = 5;

    double base = pow(q1_q2_max / q1_q2_min, 1.0 / (n_histogramas-1));
    double q1 = q1_q2_min;

    int nsteps = (int)(1e8);

    for(int i=0; i<n_histogramas; i++){
        vector <double> q = {q1, q2};
        Replicadores dos_agentes(q, lam, u, du);
        dos_agentes.transitorio(10000);
        HistoLineal histo_1(0, 1, 100);
        HistoLineal histo_2(0, 1, 100);
        for(int j=0; j<nsteps; j++){
            dos_agentes.paso();
            if(j % 1000 == 0){
                histo_1.nuevo_dato(dos_agentes.get_x()[0]);
                histo_2.nuevo_dato(dos_agentes.get_x()[1]);
            }
        }
        FILE* grafico = popen("gnuplot -persist", "w");
        histo_1.plot_densidad(grafico, true);
        histo_2.plot_densidad(grafico, true);
        q1 *= base;
    }*/
    return 0;
}
