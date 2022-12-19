#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAM 100
#define STEPS 6000
#define INTERVAL 600
#define real double

int main()
{
    FILE *saida = fopen("hipertermia_heterogeneo_np_0.txt", "w");
    if (saida == NULL)
    {
        printf("Ocorreu um erro na abertura do arquivo de saída.\n");
        return 1;
    }

    int i, j, x, y, lessX, lessY, plusX, plusY, pele, k, swap = 0, **tecido;
    char str[40];
    real pb = 1000, Ta = 37, cb = 4200, *p, *c, *Q, A = 1.3*pow(10, 6), r0 = 3.1*pow(10, -3);
    real ***u, ***u_new, hx = 0.001, ht = 0.5, *wb, **Calor_nano, *kappa, kip, kim, kjp, kjm;

    tecido = malloc (TAM*sizeof(int*));
    Q = malloc(2*sizeof(real));
    Calor_nano = malloc(TAM*sizeof(real*));
    p = malloc(2*sizeof(real));
    u = malloc(2*sizeof(real**));
    u_new = malloc(10*sizeof(real**));
    c = malloc(2*sizeof(real));
    wb = malloc(2*sizeof(real));
    kappa = malloc(2*sizeof(real));
    for (i = 0; i < TAM; i++)
    {
        tecido[i] = malloc(TAM*sizeof(int));
        Calor_nano[i] = malloc(TAM*sizeof(real));
    }
    for (i = 0; i < 10; i++)
    {
        if (i < 2)
            {
                u[i] = malloc(TAM*sizeof(real*));
                u_new[i]= malloc(TAM*sizeof(real*));
            }
        else
            u_new[i] = malloc(TAM*sizeof(real*));
    }
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < TAM; j++)
        {
            if (i < 2)
            {
                u[i][j] = malloc(TAM*sizeof(real));
                u_new[i][j] = malloc(TAM*sizeof(real));
            }
            else
                u_new[i][j] = malloc(TAM*sizeof(real));
        }
    }

    // defining constants
    p[0] = 1000;
    c[0] = 4200;
    Q[0] = 4200;
    wb[0] = 1.25*(pow(10, -3));
    kappa[0] = 0.5;
    p[1] = 1200;
    c[1] = 3600;
    Q[1] = 420;
    wb[1] = 5*(pow(10, -4));
    kappa[1] = 0.55;
    for (i = 0; i < TAM; i++)
    {
        for (j = 0; j < TAM; j++)
        {
            u[swap][i][j] = Ta;
            Calor_nano[i][j] = A*exp(-(pow((i*hx - 0.05), 2) + pow((j*hx - 0.05), 2))/pow(r0, 2));
            tecido[i][j] = ((i*hx >= 0.04 + hx/2 && i*hx <= 0.06 + hx/2) && (j*hx >= 0.04 + hx/2 && j*hx <= 0.06 + hx/2)) ? 1 : 0;
        }
    }

    // Saving the initial state
    for (i = 0; i < TAM; i++)
    {
        for (j = 0; j < TAM; j++)
            fprintf(saida, "%lf ", u[swap][i][j]);

        fprintf(saida, "\n");
    }
    fclose(saida);

{
    swap = 0;
    for (k = 0; k < STEPS; k++)
    {
      for (i = 0; i < TAM; i++)
      {
          for (j = 0; j < TAM; j++)
          {
            pele = tecido[i][j];
            x = i;
            y = j;

            lessX = (x == 0)         ? x : x - 1;
            plusX = (x == (TAM - 1)) ? x : x + 1;

            kim = (x == 0)         ?  kappa[pele] : 2*kappa[tecido[lessX][y]]*kappa[pele]/(kappa[tecido[lessX][y]] + kappa[pele]);
            kip = (x == (TAM - 1)) ?  kappa[pele] : 2*kappa[tecido[plusX][y]]*kappa[pele]/(kappa[tecido[plusX][y]] + kappa[pele]);

            lessY = (y == 0)         ?  y : y - 1;
            plusY = (y == (TAM - 1)) ?  y : y + 1;

            kjm = (y == 0)         ? kappa[pele] : 2*kappa[tecido[x][lessY]]*kappa[pele]/(kappa[tecido[x][lessY]] + kappa[pele]);
            kjp = (y == (TAM - 1)) ? kappa[pele] : 2*kappa[tecido[x][plusY]]*kappa[pele]/(kappa[pele] + kappa[tecido[x][plusY]]);

            u[!swap][x][y] = (ht/(p[pele]*c[pele]))*(
                      (kip*u[swap][plusX][j] + kjp*u[swap][x][plusY] - (kip+kim+kjp+kjm)*u[swap][i][j] + kim*u[swap][lessX][j] + kjm*u[swap][x][lessY])/(pow(hx, 2))
                      + wb[pele]*pb*cb*(Ta - u[swap][i][j]) + Q[pele] + Calor_nano[i][j]) + u[swap][i][j];
          }
      }
      swap = !swap;
      if (k % INTERVAL == 0)
        for (i = 0; i < TAM; i++)
            for (j = 0; j < TAM; j++)
                u_new[k/INTERVAL][i][j] = u[swap][i][j];
    }
}
    for (k = 0; k < 10; k++)
    {
        printf("Salvando t=%f\n", ht*k*INTERVAL);
        sprintf(str, "hipertermia_heterogeneo_np_%d.txt", k+1);
        FILE *saida = fopen(str, "w");
        for (i = 0; i < TAM; i++)
        {
            for (j = 0; j < TAM; j++)
                fprintf(saida, "%lf ", u_new[k][i][j]);

            fprintf(saida, "\n");
        }
        fclose(saida);
    }

    for (i = 0; i < TAM; i++)
    {
        free(Calor_nano[i]);
        free(tecido[i]);
    }

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < TAM; j++)
            free(u[i][j]);
        free(u[i]);
    }
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < TAM; j++)
            free(u_new[i][j]);
        free(u_new[i]);
    }
    free(tecido);
    free(Q);
    free(p);
    free(u);
    free(u_new);
    free(c);
    free(wb);
    free(Calor_nano);
    free(kappa);
    return 0;
}

