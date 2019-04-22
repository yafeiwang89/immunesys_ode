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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

// declare cell definitions here 
Cell_Definition infection_cell;
Cell_Definition scout_cell;
Cell_Definition immune_cell;

double wait_time = 1000.0;
static int scouts_left = 0 ;
double x_0, x_1, y_0, y_1 = 0.0;
double total_conc = 0.0;

static int untrained_immune_cell = 50, trained_immune_cell = 20; 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// set default cell cycle model 
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	// needed for a 2-D simulation
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0032; 

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0001; 


	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0.0002; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 1.0;
 	
 	cell_defaults.phenotype.secretion.secretion_rates[1] = 0.5; 
	cell_defaults.phenotype.secretion.uptake_rates[1] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[1] = 40.0; 

	cell_defaults.phenotype.secretion.secretion_rates[2] = 0.0; 
	cell_defaults.phenotype.secretion.uptake_rates[2] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[2] = 0.0;

	cell_defaults.phenotype.mechanics.set_relative_maximum_adhesion_distance(1.6);
	

	infection_cell = cell_defaults;
	infection_cell.type = 0;
	infection_cell.name = "invader cell";

	infection_cell.phenotype.secretion.uptake_rates[0] *= 
		parameters.doubles("infection_o2_relative_uptake"); 

	// turn on motility; 
	infection_cell.phenotype.motility.is_motile = true; 
	infection_cell.phenotype.motility.persistence_time = 5.0; 
	infection_cell.phenotype.motility.migration_speed = 0.4; 
	
	// set functions 
	infection_cell.functions.update_migration_bias = update_infection_motility;	
	infection_cell.functions.update_phenotype = update_infection_phenotype;

	infection_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 
	infection_cell.phenotype.mechanics.cell_cell_repulsion_strength = 4.0; 

	infection_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.5; 
	infection_cell.phenotype.mechanics.cell_cell_repulsion_strength = 20.0;

	create_scout_cell();
	create_immune_cell();

	return; 
}

void create_scout_cell(void){

	scout_cell = cell_defaults; 
	scout_cell.type = 1; 
	scout_cell.name = "scout cell"; 
	
	// Don't want scout cells to be adhesive.
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	scout_cell.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0001; // 0.2; 
	
	scout_cell.custom_data.add_variable("trained", "dimensionless", 0.0);
	scout_cell.custom_data.add_variable("wait_time", "dimensionless", 0.0);

	scout_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	scout_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	scout_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.saturation_densities[0] = 0.0; 


	scout_cell.phenotype.death.rates[apoptosis_model_index] = 0.0001 ; 

	// turn on motility; 
	scout_cell.phenotype.motility.is_motile = true; 
	scout_cell.phenotype.motility.persistence_time = 5.0; 
	scout_cell.phenotype.motility.migration_speed = 0.5;   

	scout_cell.phenotype.mechanics.cell_cell_adhesion_strength = 1.0; 
	scout_cell.phenotype.mechanics.cell_cell_repulsion_strength = 5.0; 

	scout_cell.functions.update_migration_bias = update_scout_motility;	
	scout_cell.functions.update_phenotype = update_scout_phenotype;
	scout_cell.functions.custom_cell_rule = scout_rule;

	return;

}

void create_immune_cell(void){

	immune_cell = cell_defaults; 
	immune_cell.type = 2; 
	immune_cell.name = "trained cell"; 
	
	// Don't want scout cells to be adhesive.
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	immune_cell.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.001; // 0.2; 
	
	immune_cell.custom_data.add_variable("wait_time", "dimensionless", 0.0);

	immune_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.saturation_densities[0] = 1.0; 

	immune_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	immune_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	immune_cell.phenotype.death.rates[apoptosis_model_index] = 0.0001 ; 

	// turn on motility; 
	immune_cell.phenotype.motility.is_motile = true; 
	immune_cell.phenotype.motility.persistence_time = 10.0; 
	immune_cell.phenotype.motility.migration_speed = 0.4;   

	immune_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 
	immune_cell.phenotype.mechanics.cell_cell_repulsion_strength = 10.0; 

	immune_cell.functions.update_migration_bias = update_immune_motility;	
	immune_cell.functions.update_phenotype = update_immune_phenotype;
	immune_cell.functions.custom_cell_rule = immune_rule;

	return;

}


void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// enable gradients calculation

	default_microenvironment_options.calculate_gradients = true;

	microenvironment.add_density( "chemoattractant", "dimensionless" ); 
	microenvironment.diffusion_coefficients[1] = 1e5; 
	microenvironment.decay_rates[1] = .1;  

	microenvironment.add_density( "mem_attractant", "dimensionless" ); 
	microenvironment.diffusion_coefficients[2] = 1e5; 
	microenvironment.decay_rates[2] = .1; 
	
	// let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = true; 

	// set Dirichlet conditions 
	//  open the boundary condition, works as a source to provide the substrate 
	default_microenvironment_options.outer_Dirichlet_conditions = true; 

	
	default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // physioxic conditions 
	default_microenvironment_options.Dirichlet_condition_vector[1] = 0; 
	default_microenvironment_options.Dirichlet_condition_vector[2] = 20; 
	
	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}

// Randomly positin the cells in the domain
void setup_tissue( void )
{

	int number_of_infection = parameters.ints("number_of_infection"); // 10;   
	int number_of_scouts = parameters.ints("number_of_scouts"); // 10;  

	std::cout << "Placing cells ... " << std::endl; 
	
	// randomly place seed cells 
	
	std::vector<double> position = {3.0, 0.0, 0.0}; 

	x_0 = default_microenvironment_options.X_range[0]; 
	x_1 = default_microenvironment_options.X_range[1];

	y_0 = default_microenvironment_options.Y_range[0]; 
	y_1 = default_microenvironment_options.Y_range[1];

	double x_range = x_1 - x_0; 
	double y_range = y_1 - y_0; 

	double relative_margin = 0.2;  
	double relative_outer_margin = 0.02;

	// Creating 10 tumor cells, 15 normal and 5 scout and assigning random position.

	Cell* pC;

	//Infection cell randomly
	for( int n=0 ; n < number_of_infection ; n++ )
	{
		position[0] = default_microenvironment_options.X_range[0] + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		position[1] = default_microenvironment_options.Y_range[0] + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 

		pC = create_cell(infection_cell); 
		pC->assign_position(position);
	}

	// Scout cells
	for( int n=0 ; n < number_of_scouts ; n++ )
	{
		position[0] = default_microenvironment_options.X_range[0] + 
				x_range*( relative_outer_margin + (1-2.0*relative_outer_margin)*UniformRandom() ); 
		
		position[1] = default_microenvironment_options.Y_range[0] + 
				y_range*( relative_outer_margin + (1-2.0*relative_outer_margin)*UniformRandom() ); 

		pC = create_cell(scout_cell);
		pC->assign_position(position); 
	}

	return; 
}


void update_infection_phenotype (Cell* pCell, Phenotype& phenotype, double dt){

		if( pCell->phenotype.death.dead == true )
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			return; 
		}

		int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );

		int o2_ind = microenvironment.find_density_index( "oxygen" );
		double o2_conc = pCell->nearest_density_vector()[o2_ind];

		int chem_ind = microenvironment.find_density_index( "chemoattractant" );

		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
		double p = pCell->state.simple_pressure;

		if (o2_conc<=1.0)
		{	
			pCell->phenotype.death.rates[apoptosis_model_index] = 0.001;
			pCell->phenotype.secretion.uptake_rates[o2_ind] = 0.00001;
		}

		else {
			// chemoattractant secretion rate
			pCell->phenotype.secretion.secretion_rates[chem_ind] = 1.0;
		}

		if (p >=0.4){
			phenotype.death.rates[apoptosis_model_index] = 0.0005;
		}

	return;

}

void update_infection_motility( Cell* pCell, Phenotype& phenotype, double dt )
{		
		if( phenotype.death.dead == true )

		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			return; 
		}

		int o2_ind = microenvironment.find_density_index( "oxygen" );
		int chem_ind = microenvironment.find_density_index( "chemoattractant" );

		double o2_conc = pCell->nearest_density_vector()[o2_ind]; 
		double chem_conc = pCell->nearest_density_vector()[chem_ind]; 

		 if (chem_conc>=20.0) {
			pCell->phenotype.motility.is_motile = false;
		}

		else if(chem_conc<10.0){
			phenotype.motility.is_motile = true;
			phenotype.motility.migration_bias = 0.5;
			phenotype.motility.migration_speed = 0.3;
			phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chem_ind);	
			normalize( &( pCell->phenotype.motility.migration_bias_direction ) );
		}

		return; 
	}

void update_scout_phenotype (Cell* pCell, Phenotype& phenotype, double dt){

	if( phenotype.death.dead == true )
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int o2_ind = microenvironment.find_density_index( "oxygen" );
	double o2_conc = pCell->nearest_density_vector()[o2_ind]; 

	if (total_conc <= 10.0)
		{
			phenotype.death.rates[apoptosis_model_index] = 0.001;
		}

		return;
	}


void update_scout_motility( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true)
			{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
			}

		int chem_index = microenvironment.find_density_index( "chemoattractant" );
		double chem_conc = pCell->nearest_density_vector()[chem_index];
		int train = pCell->custom_data.find_variable_index( "trained" );

		if ((chem_conc>0.01) && (pCell->custom_data[train] == 0.0)){
		phenotype.motility.migration_bias = 0.6;
		phenotype.motility.migration_speed = 0.6;
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chem_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );
	}

	return;
}


void scout_rule( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true)
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		Cell* pC = NULL;
		static int train = pCell->custom_data.find_variable_index( "trained" );
		static int w = pCell->custom_data.find_variable_index("wait_time");
		
		int mem_index = microenvironment.find_density_index( "mem_attractant" );
		double time = 0.0;

		std::vector<double> pos = pCell->position;

		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
		
		if(pCell->custom_data[train] == 0.0)
			{
			for( int n=0 ; n < pCell->cells_in_my_container().size() ; n++ ) {
			pC = pCell->cells_in_my_container()[n]; 

			if( (pC != pCell))
			{
			
			std::vector<double> dist = pC->position;
			dist -= pCell->position;
			double distance = norm(dist);
			int temp = pC->custom_data.find_variable_index("signal");

				if((distance <= 15.0) && (pC->type == 0))
				{	
				pCell->custom_data[train] = 1.0;
				break;
				} 

			}
		}

		}

		pCell->custom_data[w] = pCell->custom_data[w] + dt;

		if((pCell->custom_data[train] > 0.5) && (pCell->custom_data[w] >= wait_time)){
			phenotype.mechanics.cell_cell_repulsion_strength = 30.0;
			phenotype.motility.is_motile = true;
			phenotype.motility.migration_bias_direction = pCell->nearest_gradient(mem_index);	
			normalize( &( phenotype.motility.migration_bias_direction ) );
			
			phenotype.motility.migration_bias = 0.6;
			phenotype.motility.migration_speed = 0.8;
			phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.016;
			}


		if ((pos[0] >= x_1-20.0 || pos[0] <= x_0+20.0 || pos[1] >= y_1-20.0 || pos[1] <= y_0+20.0) && (pCell->custom_data[train] > 0.5)){
				
				#pragma omp critical
				scouts_left += 1;	
				printf("\n***** scouts left, %d\n", scouts_left);
				pCell->functions.update_migration_bias = NULL; 
				pCell->functions.update_phenotype = NULL;
				pCell->functions.custom_cell_rule = NULL;
				pCell->flag_for_removal();
			}

	return;
}


void update_immune_phenotype (Cell* pCell, Phenotype& phenotype, double dt){


	if( phenotype.death.dead == true )
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
		int chem_ind = microenvironment.find_density_index( "chemoattractant" );
		double chem_conc = pCell->nearest_density_vector()[chem_ind]; 
		
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

		double p = pCell->state.simple_pressure;

	if (chem_conc >= 5.0)
		{
			phenotype.death.rates[apoptosis_model_index] = 0.0001;
			phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.014;
		}

	else{
			phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0;
	}	
		return;
}

void update_immune_motility( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true)
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		int chem_index = microenvironment.find_density_index( "chemoattractant" );
		phenotype.motility.migration_bias = 0.6;
		phenotype.motility.migration_speed = 0.45;
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chem_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );


	return;
}


void immune_rule( Cell* pCell, Phenotype& phenotype, double dt ){


		if( phenotype.death.dead == true)
		{ 	
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		Cell* pC = NULL;
		std::vector<double> pos = pCell->position;

		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
		int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );

		// int o2_index = microenvironment.find_density_index( "oxygen" );

		int chem_ind = microenvironment.find_density_index( "chemoattractant" );
		double chem_conc = pCell->nearest_density_vector()[chem_ind]; 

		int mem_index = microenvironment.find_density_index( "mem_attractant" );
		
		for( int n=0 ; n < pCell->cells_in_my_container().size() ; n++ ) {

		pC = pCell->cells_in_my_container()[n]; 

		if( (pC != pCell))
		{
		
		std::vector<double> dist = pC->position;
		dist -= pCell->position;
		double distance = norm(dist);

			if(distance <= 15.0)
			{	
			if (pC->type == 0){
			pC->start_death(apoptosis_model_index);
			phenotype.death.rates[apoptosis_model_index] = 0.0005;
			//phenotype.motility.is_motile = false;
			}

			if(pCell->cells_in_my_container().size() >= 4 && pC->type == 2){
				pCell->phenotype.death.rates[apoptosis_model_index] = 9999; 
			} 

		}
	}
	
	}

	if(total_conc <= 10.0){
			phenotype.motility.migration_bias = 0.6;
			phenotype.motility.migration_speed = 0.7;
			phenotype.motility.migration_bias_direction = pCell->nearest_gradient(mem_index);	
			normalize( &( phenotype.motility.migration_bias_direction ) );
	}

	if ((pos[0] >= x_1-20.0 || pos[0] <= x_0+20.0 || pos[1] >= y_1-20.0 || pos[1] <= y_0+20.0)){
				#pragma omp critical	
				printf("\n***** scouts left, %d\n", scouts_left);
				pCell->functions.update_migration_bias = NULL; 
				pCell->functions.update_phenotype = NULL;
				pCell->functions.custom_cell_rule = NULL;
				pCell->flag_for_removal();
	}

	return;
}


// Defining my own custom coloring function
std::vector<std::string> my_coloring_function( Cell* pCell )
{
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// Live tumor cells RED, dead black
	if( pCell->phenotype.death.dead == false && pCell->type == 0)
	{
		output[1] = "red";
		output[3] = "red";
		output[0] = "red"; 
		output[2] = "red";
	}
	else if( pCell->phenotype.death.dead == true && pCell->type == 0 )
	{
		output[1] = "black";
		output[3] = "black";
		output[0] = "black"; 
		output[2] = "black";
	}

	// Scout normal cells ORANGE, dead black
	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		output[1] = "orange";
		output[3] = "orange";
		output[0] = "orange"; 
		output[2] = "orange";

	}
	else if( pCell->phenotype.death.dead == true && pCell->type == 1 )
	{
		output[1] = "black";
		output[3] = "black";
		output[0] = "black"; 
		output[2] = "black";
	}

	if( pCell->phenotype.death.dead == false && pCell->type == 2 )
	{
		output[1] = "blue";
		output[3] = "blue";
		output[0] = "blue"; 
		output[2] = "blue";

	}
	else if( pCell->phenotype.death.dead == true && pCell->type == 2 )
	{
		output[1] = "black";
		output[3] = "black";
		output[0] = "black"; 
		output[2] = "black";
	}
	
	return output; 
}


void ode_func(double dt){

	int chem_index = microenvironment.find_density_index( "chemoattractant" );
	int scout_cells, new_untrained_cells, new_trained_cells = 0;

	// collect the chemoattractant concentration from the domain (all voxels)
	#pragma omp parallel for 
	for( int k=0 ; k < microenvironment.number_of_voxels(); k++ )
	{
	total_conc += microenvironment(k)[chem_index];
	}

	printf("\nConc , %f\n", total_conc);

	if ( total_conc >= 100.0 && scouts_left > 0) {

		//spawn cells
		scout_cells = fabs( dt * (total_conc  * 0.05 * scouts_left) - (0.05 * scouts_left) - (0.02 * scouts_left * untrained_immune_cell)) ;
		new_untrained_cells =  fabs(  dt * (total_conc * 0.2 * untrained_immune_cell) - (0.03 * untrained_immune_cell) - ( 0.05 * scout_cells * untrained_immune_cell)) ;
		new_trained_cells =  fabs(dt * (total_conc * 0.3 * trained_immune_cell) - (0.034 * trained_immune_cell) + (0.003  * scouts_left * untrained_immune_cell) - (0.1 * trained_immune_cell)) ;
		
		printf("\n scout cells: %d, untrained cells: %d, trained cells %d\n", scout_cells, new_untrained_cells, new_trained_cells);
		
		if (scout_cells < 20){
			scout_cells = 10;
		}
		
		if (new_trained_cells >25 || new_trained_cells < 5){
			new_trained_cells = 15;
		}
		
		printf("\n Conc. , %fscout cells: %d, untrained cells: %d, trained cells %d\n", total_conc, scout_cells, new_untrained_cells, new_trained_cells);
		
		// call the setup tissue function again.
		untrained_immune_cell = new_untrained_cells;
		trained_immune_cell = new_trained_cells;
		scouts_left = 0;
		total_conc = 0.0;
		setup_new_cells(scout_cells, new_trained_cells);

	}
	return;
}

void setup_new_cells(int scout_cells, int trained_cells){

	int number_of_trained = trained_cells;   
	int number_of_scouts = scout_cells;  

	std::cout << "Placing new cells ... " << std::endl; 

	double xperiod = 30.0;
	double yperiod = 30.0;

	std::vector<double> position = {2.0, 1.0, 0.0}; 

	Cell* pC;

	for( int n=0 ; n < number_of_scouts ; n++ )
	{
		position[0] = default_microenvironment_options.X_range[0] + xperiod;
		position[1] = default_microenvironment_options.Y_range[0] + yperiod;
		pC = create_cell(scout_cell);
		pC->assign_position(position); 
		xperiod += xperiod;
	}

	xperiod = 30.0;
	for( int n=0 ; n < number_of_trained ; n++ )
	{
		position[0] = default_microenvironment_options.X_range[0] + xperiod;
		position[1] = default_microenvironment_options.Y_range[1] - yperiod;
		pC = create_cell(immune_cell);
		pC->assign_position(position); 
		xperiod += xperiod;
	}


	return;
}


void read_write_function(void)
{

	std::vector< std::vector<double> > cells_data;  // must a 2D matrix for writing to .mat file
	std::vector<double> template_vector3(3,0);
	cells_data.resize(1, template_vector3);
	cells_data.push_back(template_vector3);
	
	
	for( int i=0 ; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i];
		if(pCell->type == 0 && pCell->phenotype.death.dead == false)
		{
			cells_data[0][0] += 1;
		}
		else if(pCell->type == 1 && pCell->phenotype.death.dead == false)
		{
			cells_data[0][1] += 1;
		}			
		else if(pCell->type == 2 && pCell->phenotype.death.dead == false)
		{
			cells_data[0][2] += 1;
		}		
		
	}
    
	
	char filename[1024];
	sprintf( filename , "%s/new_data%08u.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
	
	// write_matlab4( cells_data , filename, "cells_num");
	write_matlab( cells_data , filename);
	
	
	return ;
	
}