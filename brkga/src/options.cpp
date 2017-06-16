/***************************************************************************
 *   copyright (C) 2005 by Marco Caserta                                   *
 *   marco dot caserta at ie dot edu                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \file options.cpp 
  \brief Read options from command line.

  Options are:
  - -f : problem instance file [default = NONE]
  - -h : help (list of all options)
*/

#include <iostream>
#include <cstdlib>

#define   TIME_LIMIT_def  1000	 //!< default wall-clock time limit
#define   ITER_LIMIT_def  500	 //!< default lagrangean iterations
#define   CORRIDOR_def    0.9
#define   NSOL_def        3
#define   PROP_ZERO_def   0.0
#define   PROP_ONE_def    0.0
#define   ADD_CUT_def     0
#define   OMEGA_def       0.0
#define   SIGMA_SQ_def    0.0

using namespace std;

extern char* _FILENAME; 	//!< name of the instance file
extern double time_limit;	//!< wall-clock time limit
extern int max_iter;		//!< max number of lagrangean iterations
extern double corridorWidthBase;
extern int nSolBase;
extern double propFixed0;
extern double propFixed1;
extern int add_oldSol_cut;
extern double Omega;
extern double sigmaSq;

/// Parse command line options
int parseOptions(int argc, char* argv[])
{
   bool setFile      = false;
   time_limit        = TIME_LIMIT_def;
   max_iter          = ITER_LIMIT_def;
   corridorWidthBase = CORRIDOR_def;
   nSolBase          = NSOL_def;
   propFixed0        = PROP_ZERO_def;
   propFixed1        = PROP_ONE_def;
   add_oldSol_cut    = ADD_CUT_def;
   Omega             = OMEGA_def;
   sigmaSq           = SIGMA_SQ_def;



   cout <<endl << "MMKP v1.0 -- MC 2015(c)" << endl;
   if (argc == 1)
   {
      cout << "No options specified. Try -h " << endl;
      return -1;
   }  
 
   int i = 0;
   while (++i < argc)
   {
      const char *option = argv[i];
      if (*option != '-')
	 return i;
      else if (*option == '\0')
	 return i;
      else if (*option == '-')
      {
	 switch (*++option)
	 {
	    case '\0':
	       return i + 1;
	    case 'f':
	       _FILENAME = argv[i+1];
	       setFile = true;
	       i++;
	       break;
	    case 't':
	       time_limit = atof(argv[i+1]);
	       i++;
	       break;
	    case 'c':
	       corridorWidthBase = atof(argv[i+1]);
	       i++;
	       break;
	    case 'n':
	       nSolBase = atoi(argv[i+1]);
	       i++;
	       break;
	    case 'z':
	       propFixed0 = atof(argv[i+1]);
	       i++;
	       break;
	    case 'u':
	       propFixed1 = atof(argv[i+1]);
	       i++;
	       break;
	    case 'a':
	       add_oldSol_cut = atoi(argv[i+1]);
	       i++;
	       break;
	    case 'i':
	       max_iter = atoi(argv[i+1]);
	       i++;
	       break;
        case 'o':
           Omega = atof(argv[i+1]);
           i++;
           break;
        case 's':
           sigmaSq = atof(argv[i+1]);
           i++;
           break;
        case 'h':
	       cout << "OPTIONS :: " << endl;
	       cout << "-f : problem instance file" << endl;
	       cout << "-t : overall time limit (real)" << endl;
	       cout << "-c : base value for corridor width" << endl;
	       cout << "-n : base nr solution sought by cplex" << endl;
	       cout << "-z : % of vars fixed to zero (hard fixed)" << endl;
	       cout << "-u : % of vars fixed to 1 (soft fixed)" << endl;
	       cout << "-a : activate cut feasible solution (1-> Yes, 0-> No)" << endl;
	       cout << "-i : max number lagrangean iterations" << endl;
	       cout << "-o : Omega (robust formulation) " << endl;
	       cout << "-s : sigma squared (item uncertainty)" << endl; 
	       cout << endl;
	       return -1;
	 }
      }
   }
 
   if (setFile)
      return 0;
   else
   {
      cout <<"Option -f is mandatory. Try ./team -h" << endl;
      return -1;
   }
}
