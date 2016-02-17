#include <stdio.h>
#include <iostream>
#include <string.h>
#include "solvers.h"
#include "program_constants.h"

int main(void)
{
    int nroots_gsl = 0;
    int rtn_gsl = 0;
    // int nroots_netlibs = 0;
    // int rtn_netlibs = 0;
    // int nroots_gems = 0;
    // int rtn_gems = 0;

    int choice = 1;
    double x_gsl[4], y_gsl[4];
    // double x_netlibs[4], y_netlibs[4];
    // double x_gems[4], y_gems[4];

    float A1 = 3.;
    float B1 = 2.;
    float H1 = 0.;
    float K1 = 0.;
    float PHI_1 = 0.;
    float A2 = 3.;
    float B2 = 1.;
    float H2 = 1.;
    float K2 = -0.5;
    float PHI_2 = 0.7853981633974483;
    float area_gsl = 0.;
    float area_netlibs = 0.;
    float area_gems = 0.;

    area_gsl = ellipse_ellipse_overlap_gsl (PHI_1, A1, B1, H1, K1, PHI_2, A2, B2, H2, K2, x_gsl, y_gsl, &nroots_gsl, &rtn_gsl, choice);
    printf("%f\n", area_gsl);
    // area_netlibs = ellipse_ellipse_overlap_netlibs (PHI_1, A1, B1, H1, K1, PHI_2, A2, B2, H2, K2, x_netlibs, y_netlibs, &nroots_netlibs, &rtn_netlibs);
    // printf("%f\n", area_netlibs);
    // area_gems = ellipse_ellipse_overlap_gems (PHI_1, A1, B1, H1, K1, PHI_2, A2, B2, H2, K2, x_gems, y_gems, &nroots_gems, &rtn_gems);
    // printf("%f\n", area_gems);
    
}
