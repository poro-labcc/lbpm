/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * Multi-relaxation time LBM Model
 */
#include "analysis/rev_analysis.h"
#include "models/MRTModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"
#include <cmath>
#include <algorithm>

void REVfunc::PoreSize(double &average_pore_size, ScaLBL_MRTModel &MRT) {
    int iter = 1;
    int Nx = MRT.Nx, Ny = MRT.Ny, Nz = MRT.Nz;

    std::vector<std::vector<std::vector<double>>> distance_updated(
        Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    for (int i = 0; i < Nx; i++) //copy Distance values to distance_updated
    {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                distance_updated[i][j][k] = MRT.Distance(i, j, k);
            }
        }
    }

    bool WriteHeader = false;
    FILE *log_file = fopen("pore_size.csv", "r");
    if (log_file != NULL)
        fclose(log_file);
    else
        WriteHeader = true;

    if (WriteHeader) {
        log_file = fopen("pore_size.csv", "a+");
        fprintf(log_file, "iteration pore_size\n");
        fclose(log_file);
    }

    double current_max = 0, prev_max = 1000000;

    while (prev_max > 0) {

        current_max = 0;

        // find current maximum
        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    double value = MRT.Distance(i, j, k);
                    if (value > current_max && value < prev_max) {
                        current_max = value;
                    }
                }
            }
        }

        if (current_max > 0)
            break;

        // after finding it, proceeds to change values inside sphere
        for (int k = 1; k < Nz - 1; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    double value = MRT.Distance(i, j, k);
                    if (value == current_max) {
                        double start_x =
                            std::max(1, static_cast<int>(i - current_max - 1));
                        double end_x = std::min(
                            Nx - 1, static_cast<int>(i + current_max + 2));

                        double start_y =
                            std::max(1, static_cast<int>(j - current_max - 1));
                        double end_y = std::min(
                            Ny - 1, static_cast<int>(j + current_max + 2));

                        double start_z =
                            std::max(1, static_cast<int>(k - current_max - 1));
                        double end_z = std::min(
                            Nz - 1,
                            static_cast<int>(
                                k + current_max +
                                2)); // o primeiro +1 para arredondar para cima, o segundo pois o for é < e não <=

                        for (int n = start_z; n < end_z; n++) {
                            for (int m = start_y; m < end_y; m++) {
                                for (int l = start_x; l < end_x; l++) {
                                    double value2 = distance_updated[l][m][n];
                                    if (value2 > 0 && value2 < current_max &&
                                        (i - l) * (i - l) + (j - m) * (j - m) +
                                                (k - n) * (k - n) <=
                                            current_max * current_max)
                                        distance_updated[l][m][n] = current_max;
                                }
                            }
                        }
                    }
                }
            }
        }

        log_file = fopen("pore_size.csv", "a+");
        fprintf(log_file, "%d %f\n", iter, 2 * current_max);
        fclose(log_file);

        prev_max = current_max;
        iter++;
    }

    // Calculate the average pore size
    double sum = 0.0;
    int count = 0;
    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                if (distance_updated[i][j][k] > 0) {
                    sum += distance_updated[i][j][k];
                    count++;
                }
            }
        }
    }
    average_pore_size = 2 * sum / count; //convert radius to diameter

    log_file = fopen("pore_size.csv", "a+");
    fprintf(log_file, "\n\n average_pore_size %f", average_pore_size);
    fclose(log_file);

    std::cout << "\n\n\nAverage Pore Size: " << average_pore_size << "\n\n\n"
              << std::endl;
}

void REVfunc::RevAnalysis(ScaLBL_MRTModel &MRT) {
    double average_pore_size = 1.0;

    PoreSize(average_pore_size, MRT);

    int Nx = MRT.Nx, Ny = MRT.Ny, Nz = MRT.Nz;
    double rev_x_perm = -1.0, rev_x_poro = -1.0, rev_x_tort = -1.0;
    int iter_step = 0, side_step = static_cast<int>(average_pore_size);
    double gamma = 0.1;
    double h = MRT.Dm->voxel_length;
    double mu = (MRT.tau - 0.5) / 3.f;
    int samples = 20; //samples to be used in deterministic standard deviation
    int stat_samples = 100; //samples to be used in statistical approach
    int x_side = 2 * static_cast<int>(average_pore_size); //starting x-side

    // Arrays for porosity, permeability and subdomain size

    int vector_size_x = static_cast<int>(std::ceil(MRT.Nx / side_step)) + 1;
    std::vector<double> sub_x_size(vector_size_x);
    std::vector<double> perm(vector_size_x);
    std::vector<double> poro(vector_size_x);
    std::vector<double> tort(vector_size_x);

    while (x_side < Nx - 2) {
        bool WriteHeader = false;
        FILE *log_file = fopen("rev_analysis.csv", "r");
        if (log_file != NULL)
            fclose(log_file);
        else
            WriteHeader = true;

        if (WriteHeader) {
            log_file = fopen("rev_analysis.csv", "a+");
            fprintf(
                log_file, "iter_step average_pore_size x_side "
                          "porosity rev_poro RE_poro CC_poro stat_poro stat_poro_mean stat_poro_max stat_poro_min "
                          "permeability rev_perm RE_perm CC_perm stat_perm stat_perm_mean stat_perm_max stat_perm_min "
                          "tortuosity rev_tort RE_tort CC_tort stat_tort stat_tort_mean stat_tort_max stat_tort_min\n");
            fclose(log_file);
        }

        // Find deterministic subdomain dimensions
        int y_side = (Ny - 2) * x_side / (Nx - 2);
        int z_side = (Nz - 2) * x_side / (Nx - 2);

        int start_x = (Nx - 2 - x_side) / 2 + 1;
        int start_y = (Ny - 2 - y_side) / 2 + 1;
        int start_z = (Nz - 2 - z_side) / 2 + 1;

        int end_x = (Nx - 2 + x_side) / 2 + 1;
        int end_y = (Ny - 2 + y_side) / 2 + 1;
        int end_z = (Nz - 2 + z_side) / 2 + 1;

        start_x = std::max(1, start_x);
        end_x = std::min(Nx-2, end_x);
        start_y = std::max(1, start_y);
        end_y = std::min(Ny-2, end_y);
        start_z = std::max(1, start_z);
        end_z = std::min(Nz-2, end_z);


        double count = 0.0;
        double vax = 0.0, vay = 0.0, vaz = 0.0;
        double v_tort=0.0;

        for (int k = start_z; k <= end_z; k++) { 
            for (int j = start_y; j <= end_y; j++) {
                for (int i = start_x; i <= end_x; i++) {
                    if (MRT.Distance(i, j, k) > 0) {
                        vax += MRT.Velocity_x(i, j, k);
                        vay += MRT.Velocity_y(i, j, k);
                        vaz += MRT.Velocity_z(i, j, k);
                        v_tort += sqrt(MRT.Velocity_x(i, j, k)*MRT.Velocity_x(i, j, k) + MRT.Velocity_y(i, j, k)*MRT.Velocity_y(i, j, k) + MRT.Velocity_z(i, j, k)*MRT.Velocity_z(i, j, k));
                        count += 1.0;
                    }
                }
            }
        }

        if (count == 0)
        {
            vax = 0.0;
            vay = 0.0;
            vaz = 0.0;
        } else {
        vax /= count;
        vay /= count;
        vaz /= count;}

        // Default to z direction
        double force_mag = sqrt(MRT.Fx * MRT.Fx + MRT.Fy * MRT.Fy + MRT.Fz * MRT.Fz);
        double dir_x = MRT.Fx / force_mag;
        double dir_y = MRT.Fy / force_mag;
        double dir_z = MRT.Fz / force_mag;
        if (force_mag == 0.0) {
            // default to z direction
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
            force_mag = 1.0;
        }

        double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);

        double current_poro = count / ((end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1)); //porosity
        double current_perm = 1013 * h * h * mu * current_poro * flow_rate / force_mag; //permeability
        double current_tort = 1e100;

        if (vaz != 0)
            current_tort = v_tort / (vaz*count) - 1;

        // Saving data
        sub_x_size[iter_step] = x_side;
        poro[iter_step] = current_poro;
        perm[iter_step] = current_perm;
        tort[iter_step] = current_tort;
        


        // Porosity


        
        //Porosity deterministic approach

        // Convergence criteria
        double CC = 10.0, RE = 10.0;
        double mean = 0.0;
        double strd_dev = 0.0;


        if (iter_step >= 1 && iter_step < samples - 1)
        {
            for (int l = 0; l <= iter_step; l++)
                mean += poro[l] / (iter_step + 1);

            for (int l = 0; l <= iter_step; l++)
                strd_dev += (poro[l] - mean) * (poro[l] - mean);

            strd_dev = sqrt(strd_dev / (iter_step + 1));

            if(mean != 0)
            CC = strd_dev / mean;

            if((poro[iter_step] + poro[(iter_step - 1)]) != 0)
            RE = 2 * fabs((poro[iter_step] - poro[(iter_step - 1)]) /
                          (poro[iter_step] + poro[(iter_step - 1)]));
        }

        if (iter_step >= samples - 1) {

            for (int l = iter_step + 1 - samples; l <= iter_step; l++)
                mean += poro[l] / samples;

            for (int l = iter_step + 1 - samples; l <= iter_step; l++)
                strd_dev += (poro[l] - mean) * (poro[l] - mean);

            strd_dev = sqrt(strd_dev / samples);

            if(mean != 0)
            CC = strd_dev / mean;

            if((poro[iter_step] + poro[(iter_step - 1)]) != 0)
            RE = 2 * fabs((poro[iter_step] - poro[(iter_step - 1)]) /
                          (poro[iter_step] + poro[(iter_step - 1)]));
        }

        if (RE < gamma && CC < gamma && rev_x_poro == -1) {
            rev_x_poro = x_side * h;
        }


        //Porosity statistical approach

        double stat_CC = 10;
        double stat_min = 1e100;
        double stat_max = -1e100;
        mean = 0.0;
        strd_dev = 0.0;

        if (rev_x_poro != -1) {

            std::vector<double> stat_poro(stat_samples);

            // Ensure sub_x_size[iter_step] is a positive number
            int sub_x = static_cast<int>(sub_x_size[iter_step]);
            int sub_y = std::floor((sub_x) * (Ny - 2) / (Nx - 2)) + 1;
            int sub_z = std::floor((sub_x) * (Nz - 2) / (Nx - 2)) + 1;

            printf("\n\nNy = %d, sub_y = %d\n\n", Ny, sub_y);

            for (int l = 0; l < stat_samples; l++) {
                double stat_current_pore_count = 0.0;

                // Calculate the range for the starting positions
                int range_x = (Nx - 2) - sub_x +
                              1; // +1 to include the last valid start position
                int range_y = (Ny - 2) - sub_y + 1;
                int range_z = (Nz - 2) - sub_z + 1;

                // Randomly select start positions within valid ranges
                int start_x = 1 + rand() % range_x;
                int start_y = 1 + rand() % range_y;
                int start_z = 1 + rand() % range_z;

                // Calculate end positions
                int end_x = start_x + sub_x - 1;
                int end_y = start_y + sub_y - 1;
                int end_z = start_z + sub_z - 1;

                start_x = std::max(1, start_x);
                end_x = std::min(Nx-2, end_x);
                start_y = std::max(1, start_y);
                end_y = std::min(Ny-2, end_y);
                start_z = std::max(1, start_z);
                end_z = std::min(Nz-2, end_z);


                // Count pores
                for (int k = start_z; k <= end_z; k++) {
                    for (int j = start_y; j <= end_y; j++) {
                        for (int i = start_x; i <= end_x; i++) {
                            if (MRT.Distance(i, j, k) > 0) {
                                stat_current_pore_count++;
                            }
                        }
                    }
                }

                // Compute porosity
                    stat_poro[l] = stat_current_pore_count / ((end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1));
            }

            for (int l = 0; l < stat_samples; l++) {
                mean += stat_poro[l] / (stat_samples);
                if (stat_poro[l] > stat_max)
                    stat_max = stat_poro[l];
                if (stat_poro[l] < stat_min)
                    stat_min = stat_poro[l];
            }

            for (int l = 0; l < stat_samples; l++)
                strd_dev += (stat_poro[l] - mean) * (stat_poro[l] - mean);

            strd_dev = sqrt(strd_dev / (stat_samples));

            if(mean != 0)
            stat_CC = strd_dev / mean;
        }

        log_file = fopen("rev_analysis.csv", "a");
        fprintf(log_file, "%d %f %f %f %f %f %f %f %f %f %f ", iter_step, average_pore_size * h, sub_x_size[iter_step] * h,
         poro[iter_step], rev_x_poro, RE, CC, stat_CC, mean, stat_max, stat_min);
        fclose(log_file);
        // }




        // Permeability


        //Permeability deterministic approach

        CC = 10;
        RE = 10;
        mean = 0.0;
        strd_dev = 0.0;

        if (iter_step >= 1 && iter_step < (samples - 1))
        {

            for (int l = 0; l < iter_step + 1; l++)
                mean += perm[l] / (iter_step + 1);

            for (int l = 0; l <= iter_step; l++)
                strd_dev += (perm[l] - mean) * (perm[l] - mean);

            strd_dev = sqrt(strd_dev / (iter_step + 1));

            if(mean !=0)
            CC = strd_dev / mean;

            if((perm[iter_step] + perm[(iter_step - 1)]) != 0)
            RE = 2 * fabs((perm[iter_step] - perm[(iter_step - 1)]) / (perm[iter_step] + perm[(iter_step - 1)]));
        }

        if (iter_step >= samples - 1)
        {

            for (int l = iter_step + 1 - samples; l <= iter_step; l++)
                mean += perm[l] / samples;

            for (int l = iter_step + 1 - samples; l <= iter_step; l++)
                strd_dev += (perm[l] - mean) * (perm[l] - mean);

            strd_dev = sqrt(strd_dev / samples);

            if(mean !=0)
            CC = strd_dev / mean;

            if((perm[iter_step] + perm[(iter_step - 1)]) != 0)
            RE = 2 * fabs((perm[iter_step] - perm[(iter_step - 1)]) / (perm[iter_step] + perm[(iter_step - 1)]));
        }

        if (RE < gamma && CC < gamma && rev_x_perm == -1) {
            rev_x_perm = x_side * h;
        }


        // Permeability statistical approach

        stat_CC = -1.0;
        stat_min = 1e100;
        stat_max = -1e100;
        mean = 0.0;
        strd_dev = 0.0;

        if (rev_x_perm != -1) {
            std::vector<double> stat_perm(stat_samples);

            // Ensure sub_x_size[iter_step] is a positive number
            int sub_x = static_cast<int>(sub_x_size[iter_step]);
            int sub_y = std::floor((sub_x) * (Ny - 2) / (Nx - 2)) + 1;
            int sub_z = std::floor((sub_x) * (Nz - 2) / (Nx - 2)) + 1;

            printf("\n\nNy = %d, sub_y = %d\n\n", Ny, sub_y);

            for (int l = 0; l < stat_samples; l++) {

                // Calculate the range for the starting positions
                int range_x = (Nx - 2) - sub_x +
                              1; // +1 to include the last valid start position
                int range_y = (Ny - 2) - sub_y + 1;
                int range_z = (Nz - 2) - sub_z + 1;

                // Randomly select start positions within valid ranges
                int start_x = 1 + rand() % range_x;
                int start_y = 1 + rand() % range_y;
                int start_z = 1 + rand() % range_z;

                // Calculate end positions
                int end_x = start_x + sub_x - 1;
                int end_y = start_y + sub_y - 1;
                int end_z = start_z + sub_z - 1;


                start_x = std::max(1, start_x);
                end_x = std::min(Nx-2, end_x);
                start_y = std::max(1, start_y);
                end_y = std::min(Ny-2, end_y);
                start_z = std::max(1, start_z);
                end_z = std::min(Nz-2, end_z);


                vax = 0.0;
                vay = 0.0;
                vaz = 0.0;
                count = 0.0;

                for (int k = start_z; k <= end_z; k++) {
                    for (int j = start_y; j <= end_y; j++) {
                        for (int i = start_x; i <= end_x; i++) {
                            if (MRT.Distance(i, j, k) > 0) {
                                vax += MRT.Velocity_x(i, j, k);
                                vay += MRT.Velocity_y(i, j, k);
                                vaz += MRT.Velocity_z(i, j, k);
                                count += 1.0;
                            }
                        }
                    }
                }

                if (count == 0)
                {
                    vax = 0.0;
                    vay = 0.0;
                    vaz = 0.0;
                } else {
                    vax /= count;
                    vay /= count;
                    vaz /= count;}

                // Default to z direction
                double force_mag = sqrt(MRT.Fx * MRT.Fx + MRT.Fy * MRT.Fy + MRT.Fz * MRT.Fz);
                double dir_x = MRT.Fx / force_mag;
                double dir_y = MRT.Fy / force_mag;
                double dir_z = MRT.Fz / force_mag;
                if (force_mag == 0.0) {
                    // default to z direction
                    dir_x = 0.0;
                    dir_y = 0.0;
                    dir_z = 1.0;
                    force_mag = 1.0;
                }

                double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);

                current_poro = count / ((end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1));
                stat_perm[l] = 1013 * h * h * mu * current_poro * flow_rate / force_mag;
            }

            for (int l = 0; l < stat_samples; l++) {
                mean += stat_perm[l] / (stat_samples);
                if (stat_perm[l] > stat_max)
                    stat_max = stat_perm[l];
                if (stat_perm[l] < stat_min)
                    stat_min = stat_perm[l];
            }

            for (int l = 0; l < stat_samples; l++)
                strd_dev += (stat_perm[l] - mean) * (stat_perm[l] - mean);

            strd_dev = sqrt(strd_dev / (stat_samples));

            stat_CC = strd_dev / mean;
        }

        log_file = fopen("rev_analysis.csv", "a");
        fprintf(log_file, "%f %f %f %f %f %f %f %f ", perm[iter_step], rev_x_perm, RE, CC, stat_CC, mean, stat_max, stat_min);
        fclose(log_file);

        // tortuosity

        RE = 10;
        CC = 10;

        if (iter_step >= 1 && iter_step <= samples - 1)
        {
            mean = 0;
            strd_dev = 0;
            double tort_count = 0;

            for (int l = 0; l <= iter_step; l++)
                if(tort[l] != 1e100){
                mean += tort[l];
                tort_count++;}

            if(tort_count != 0)
            mean /= tort_count;

            for (int l = 0; l <= iter_step; l++)
                if(tort[l] != 1e100)
                strd_dev += (tort[l] - mean) * (tort[l] - mean);

            strd_dev = sqrt(strd_dev / tort_count);

            if(mean != 0)
            CC = strd_dev / mean;
            if((tort[iter_step] + tort[(iter_step - 1)]) != 0)
            RE = 2 * fabs((tort[iter_step] - tort[(iter_step - 1)]) / (tort[iter_step] + tort[(iter_step - 1)]));
        }

        if (iter_step > samples - 1)
        {
            mean = 0;
            strd_dev = 0;
            double tort_count = 0;

            for (int l = iter_step + 1 - samples; l <= iter_step; l++)
                if(tort[l] != 1e100){
                mean += tort[l] / samples;
                tort_count++;}

            for (int l = iter_step + 1 - samples; l <= iter_step; l++)
                if(tort[l]!= 1e100)
                strd_dev += (tort[l] - mean) * (tort[l] - mean);

            strd_dev = sqrt(strd_dev / tort_count);

            if (mean != 0)
            CC = strd_dev / mean;

            if((tort[iter_step] + tort[(iter_step - 1)]) != 0)
            RE = 2 * fabs((tort[iter_step] - tort[(iter_step - 1)]) / (tort[iter_step] + tort[(iter_step - 1)]));
        }

        if (RE < gamma && CC < gamma && rev_x_tort == -1)
        {
            rev_x_tort = x_side * h;
        }

        stat_CC = -1.0;
        stat_min = 1e100;
        stat_max = -1e100;
        mean = 0.0;
        strd_dev = 0.0;

        // statistical tortuosity

        if (rev_x_tort != -1)
        {

            std::vector<double> stat_tort(stat_samples);

            // Ensure sub_x_size[iter_step] is a positive number
            int sub_x = static_cast<int>(sub_x_size[iter_step]);
            int sub_y = std::floor((sub_x) * (Ny - 2) / (Nx - 2)) + 1;
            int sub_z = std::floor((sub_x) * (Nz - 2) / (Nx - 2)) + 1;

            printf("\n\nNy = %d, sub_y = %d\n\n", Ny, sub_y);

            for (int l = 0; l < stat_samples; l++) {

                // Calculate the range for the starting positions
                int range_x = (Nx - 2) - sub_x +
                              1; // +1 to include the last valid start position
                int range_y = (Ny - 2) - sub_y + 1;
                int range_z = (Nz - 2) - sub_z + 1;

                // Randomly select start positions within valid ranges
                int start_x = 1 + rand() % range_x;
                int start_y = 1 + rand() % range_y;
                int start_z = 1 + rand() % range_z;

                // Calculate end positions
                int end_x = start_x + sub_x - 1;
                int end_y = start_y + sub_y - 1;
                int end_z = start_z + sub_z - 1;

                start_x = std::max(1, start_x);
                end_x = std::min(Nx-2, end_x);
                start_y = std::max(1, start_y);
                end_y = std::min(Ny-2, end_y);
                start_z = std::max(1, start_z);
                end_z = std::min(Nz-2, end_z);


                v_tort = 0.0;
                vaz = 0.0;
               

                for (int k = start_z; k <= end_z; k++) // LEMBRANDO: camadas extra
                {
                    for (int j = start_y; j <= end_y; j++)
                    {
                        for (int i = start_x; i <= end_x; i++)
                        {
                            if (MRT.Distance(i, j, k) > 0)
                            {
                                v_tort += sqrt(MRT.Velocity_x(i, j, k)*MRT.Velocity_x(i, j, k) + MRT.Velocity_y(i, j, k)*MRT.Velocity_y(i, j, k) + MRT.Velocity_z(i, j, k)*MRT.Velocity_z(i, j, k));
                                vaz += MRT.Velocity_z(i,j,k);
                            }
                        }
                    }
                }

                if (vaz > 0)
                    stat_tort[l]= v_tort / vaz - 1;
                else
                    stat_tort[l] = 1e100;
            }

            double stat_tort_count = 0;

            for (int l = 0; l < stat_samples; l++)
            {
                if (stat_tort[l] != 1e100){
                mean += stat_tort[l];
                stat_tort_count++;
                if (stat_tort[l] > stat_max)
                    stat_max = stat_tort[l];
                if (stat_tort[l] < stat_min)
                    stat_min = stat_tort[l];
                }
            }
            mean = mean / stat_tort_count;

            for (int l = 0; l < stat_samples; l++)
                if (stat_tort[l] != 1e100)
                    strd_dev += (stat_tort[l] - mean) * (stat_tort[l] - mean);

            strd_dev = sqrt(strd_dev / (stat_tort_count));

            stat_CC = strd_dev / mean;
            mean +=1;
        }

        // if (rank == 0)
        // { // aqui, nada mais está em voxels
            log_file = fopen("rev_analysis.csv", "a");
            fprintf(log_file, "%f %f %f %f %f %f %f %f\n", tort[iter_step] + 1, rev_x_tort, RE, CC, stat_CC, mean, stat_max+1, stat_min+1);
            fclose(log_file);
        //}

        printf(" %d\n", iter_step);
        x_side += side_step;
        iter_step++;
    }


        double vax = 0.0, vay = 0.0, vaz = 0.0;
        double count = 0.0;
        double v_tort = 0.0;

        for (int k = 1; k <= Nx-2; k++) { 
            for (int j = 1; j <= Ny-2; j++) {
                for (int i = 1; i <= Nz-2; i++) {
                    if (MRT.Distance(i, j, k) > 0) {
                        vax += MRT.Velocity_x(i, j, k);
                        vay += MRT.Velocity_y(i, j, k);
                        vaz += MRT.Velocity_z(i, j, k);
                        v_tort += sqrt(MRT.Velocity_x(i, j, k)*MRT.Velocity_x(i, j, k) + MRT.Velocity_y(i, j, k)*MRT.Velocity_y(i, j, k) + MRT.Velocity_z(i, j, k)*MRT.Velocity_z(i, j, k));
                        count += 1.0;
                    }
                }
            }
        }

        if (count == 0)
        {
            vax = 0.0;
            vay = 0.0;
            vaz = 0.0;
        } else {
        vax /= count;
        vay /= count;
        vaz /= count;}

        // Default to z direction
        double force_mag = sqrt(MRT.Fx * MRT.Fx + MRT.Fy * MRT.Fy + MRT.Fz * MRT.Fz);
        double dir_x = MRT.Fx / force_mag;
        double dir_y = MRT.Fy / force_mag;
        double dir_z = MRT.Fz / force_mag;
        if (force_mag == 0.0) {
            // default to z direction
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
            force_mag = 1.0;
        }

        double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);

        double current_poro = count / ((Nx-2) * (Ny-2) * (Nz-2)); //porosity
        double current_perm = 1013 * h * h * mu * current_poro * flow_rate / force_mag; //permeability
        double current_tort = 1e100;
                if (vaz != 0)
            current_tort = v_tort / (vaz*count);


            FILE *log_file = fopen("rev_analysis.csv", "a");
            fprintf(log_file, "- %f %d %f %f - - - - - - %f %f - - - - - - %f %f - - - - - -", average_pore_size, Nx-2, current_poro, rev_x_poro, current_perm, rev_x_perm, current_tort, rev_x_tort);
            fclose(log_file);

}
