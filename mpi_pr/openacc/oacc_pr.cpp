#include <math.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "oacc_pr.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <typeinfo>

using namespace std;

void get_distances(double **x, double **y, double **z, const int nframes, const int natoms, std::vector<std::vector<int> >& hist, const int nbins, const double bin_width) {

    int i,j,k,l ;
    unsigned long long npairs ;
    double x1, y1, z1, x2, y2, z2, sdist ;
    double dx2, dy2, dz2 ;
    double this_low_bin, this_high_bin ;
    unsigned long long local_count ; 
    unsigned long long local_hist[nbins] ; 
    unsigned long long m ;

    std::ostringstream sstream;
    std::string remark = "#oacc_pr_output";
    const std::string filename = "dum_oacc.txt";
    
    std::ofstream outfile(filename.c_str()) ;
    outfile << remark << std::endl;
  
    printf("oacc: %d\n", nframes) ;
    printf("oacc: %d\n", natoms) ;
    printf("oacc: %d\n", natoms) ;
    printf("oacc: frame 0 : atom 0\n") ;
    printf("oacc: x[0][0] = %f\n", x[0][0]) ;
    printf("oacc: y[0][0] = %f\n", y[0][0]) ;
    printf("oacc: z[0][0] = %f\n", z[0][0]) ;
    printf("oacc: frame 0 : atom 1\n") ;
    printf("oacc: x[0][1] = %f\n", x[0][1]) ;
    printf("oacc: y[0][1] = %f\n", y[0][1]) ;
    printf("oacc: z[0][1] = %f\n", z[0][1]) ;
    printf("oacc: frame 1 : atom 0\n") ;
    printf("oacc: x[1][0] = %f\n", x[1][0]) ;
    printf("oacc: y[1][0] = %f\n", y[1][0]) ;
    printf("oacc: z[1][0] = %f\n", z[1][0]) ;
    printf("oacc: frame 1 : atom 1\n") ;
    printf("oacc: x[1][1] = %f\n", x[1][1]) ;
    printf("oacc: y[1][1] = %f\n", y[1][1]) ;
    printf("oacc: z[1][1] = %f\n", z[1][1]) ;

    npairs = (natoms * (natoms - 1))/2 ;
    double dist[npairs] ;
    double local_dist[npairs] ; 
    
    std::cout << "nframes * npairs = " << nframes*npairs << std::endl ;
    std::cout << "nf * np = " << nframes * npairs << std::endl ;
    std::cout << "nf = " << nframes << std::endl ;
    std::cout << "np = " << npairs << std::endl ;

    //output_file.open("dum.txt") ;

    std::cout << "starting parallel loops" << std::endl ; 
    #pragma acc data copyin(x[nframes][natoms], y[nframes][natoms], z[nframes][natoms]) copy(dist[npairs], local_dist[npairs], local_hist[nbins])
    //#pragma acc data copyin(coor[nframes][natoms][3]) copy(dist[npairs], local_dist[npairs], local_hist[nbins])
    {
    for(i=0 ; i < nframes ; i++){
        for(m=0; m < npairs ; m++) { local_dist[m] = 0.0 ; }
        for(j=0; j < nbins ; j++) { local_hist[j] = 0 ; }
        std::cout << "." << std::flush ;
        #pragma acc parallel loop
        {
        for(j=0 ; j < natoms-1 ; j++){
            x1 = x[i][j] ;
            y1 = y[i][j] ;
            z1 = z[i][j] ;
            #pragma acc loop 
            {
            for(k=j+1 ; k < natoms ; k++){
                x2 = x[i][k] ;
                y2 = y[i][k] ;
                z2 = z[i][k] ;
                dx2 = (x1 - x2) * (x1 - x2) ;
                dy2 = (y1 - y2) * (y1 - y2) ;
                dz2 = (z1 - z2) * (z1 - z2) ;
                sdist = sqrt(dx2 + dy2 + dz2) ;
                local_count = ((j*natoms)-((j*(j+1))/2))+k-(j+1) ;
                dist[local_count] += sdist ; 
                local_dist[local_count] = sdist ; 

                #pragma acc loop 
                for(l=0; l<nbins; l++){
                    this_low_bin = double(l)*bin_width ;
                    this_high_bin = this_low_bin + bin_width ;
                    if(sdist > this_low_bin && sdist <= this_high_bin){
                        local_hist[l] += 1 ;
                        break;
                    }
                } // end of l-loop

            } // end of k-loop
            } // end of pragma acc loop
        } // end of j-loop
        } // end of pragma acc parallel loop
        
        if(i<nframes){ 
            for(k=0 ; k < nbins ; k++){
                hist[i][k] = local_hist[k] ;
            }
        }

    } // end of i-loop

    } // end of pragma acc data 
    
    for(i=0 ; i < nframes ; i++){
        for(k=0 ; k < nbins ; k++){
            outfile << k << "\t" << hist[i][k] << std::endl;
        }
        outfile << std::endl ;
    }

    outfile.close();

    std::cout << std::endl << "leaving oacc" << std::endl ;

    return ;
}
