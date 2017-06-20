/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 *      Adapted : marco caserta (for the mmkp problem)
 */



#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <fstream>

#include <sstream>
#include "SampleDecoder.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>

SampleDecoder::SampleDecoder()  { }
SampleDecoder::~SampleDecoder() { }

using namespace std;


/// FUNCTIONS DEFINITION ==========================
/// END FUNCTIONS DEFINITION ==========================

double select_corridor_width_base(double r)
{
    if (r < 0.20)
        return 0.8;
    else if(r < 0.40)
        return 0.85;
    else if (r < 0.60)
        return 0.90;
    else if (r < 0.8)
        return 0.95;
    else
        return 0.98;
}


double select_zero_base(double r)
{
    if (r < 0.25)
        return 0.6;
    else if(r < 0.50)
        return 0.7;
    else if (r < 0.75)
        return 0.8;
    else
        return 0.9;
}

double select_one_base(double r)
{
    if (r < 0.20)
        return 0.05;
    else if(r < 0.40)
        return 0.1;
    else if (r < 0.60)
        return 0.15;
    else if (r < 0.8)
        return 0.20;
    else
        return 0.25;
}

double select_nr_columns_DW(double r)
{
    if (r < 0.20)
        return 50;
        else if (r < 0.40)
            return 100;
        else if (r < 0.60)
            return 200;
        else if (r < 0.80)
            return 400;
        else
            return 1000;
}


double SampleDecoder::decode(const std::vector< double >& chromosome) const 
{
    // decoding (chromosome of length "n"):
    double corridorWidthBase = select_corridor_width_base(chromosome[0]);
    double propFixed0        = select_zero_base(chromosome[1]);
    double propFixed1        = select_one_base(chromosome[2]);
    int    nSolInPool        = select_nr_columns_DW(chromosome[3]);


    string sBase = "python ../clspBenders.py -i ../../benders/data/tr24-15 -a 3";
    stringstream s1;
    s1 << corridorWidthBase;
    stringstream s2;
    s2 << propFixed0;
    stringstream s3;
    s3 << propFixed1;
    stringstream s4;
    s4 << nSolInPool;

    string commLine =  sBase + " -c " + s1.str() + " -z " + s2.str() + " -o " 
        + s3.str() + " -p " + s4.str();

    const char * cLine = commLine.c_str();
    double outCL = system(cLine);

    string temp;
    double score;

    ifstream fsol("solution.txt", ios::in);
    fsol >> temp >> temp >> temp >> score;
    fsol.close();


    // REM. : It is minimizing...
    cout << "################ score of this chromosome is " << score << endl;

    double result = system("cat solution.txt >> summary.txt");

    // return -score; // choose this if maximazing
    return score; // choose this if minimizing
}

