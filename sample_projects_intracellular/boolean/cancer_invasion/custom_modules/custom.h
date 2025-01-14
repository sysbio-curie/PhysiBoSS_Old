/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;
#include "custom_main.h"

// setup functions to help us along 

void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// custom pathology coloring function 

std::vector<std::string> my_coloring_function( Cell* );

// custom functions can go here 

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius);
std::vector<std::vector<double>> create_cell_disc_positions(double cell_radius, double disc_radius);

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );

/** \brief Get the current value of integrin strength */
inline double get_integrin_strength( double percent )
{ return current_value( PhysiCell::parameters.doubles("ecm_adhesion_min"), PhysiCell::parameters.doubles("ecm_adhesion_max"), percent ); };

/** \brief Get the current value of motility coefficient */
inline double get_motility_amplitude( double percent )
{ return current_value(PhysiCell::parameters.doubles("motility_amplitude_min"), PhysiCell::parameters.doubles("motility_amplitude_max"), percent ); };

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ); 

/** \brief Return amount of contact with ECM */
inline double contact_ecm(Cell* pCell)
{ double ecm_contact = pCell->custom_data["ecm_contact"]; return ecm_contact / pCell->phenotype.geometry.radius ; };

void add_ecm_interaction( Cell* pCell, int index_ecm, int index_voxel );
void pre_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt);
void post_update_intracellular(Cell* pCell, Phenotype& phenotype, double dt);

std::vector<std::string> my_coloring_function( Cell* );
std::vector<std::string> ECM_coloring_function( Cell*);
std::vector<std::string> phase_coloring_function( Cell* );
std::vector<std::string> node_coloring_function( Cell* );
std::vector<std::string> my_coloring_function_for_stroma( double concentration, double max_conc, double min_conc );
void color_node(Cell* pCell);