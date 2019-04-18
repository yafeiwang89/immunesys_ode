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

float low = 0.0 ;
float high = 0.0 ;
double wait_time = 1000.0;
static int scouts_left = 0 ;
double x_0, x_1, y_0, y_1 = 0.0;

int untrained_immune_cell = 100, trained_immune_cell = 20; 

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

	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.01; 

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 


	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0.0001; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 1000.0;
 	
 	cell_defaults.phenotype.secretion.secretion_rates[1] = 10.0; 
	cell_defaults.phenotype.secretion.uptake_rates[1] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[1] = 50.0; 

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
	//infection_cell.functions.update_phenotype = update_infection_phenotype;

	infection_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 
	infection_cell.phenotype.mechanics.cell_cell_repulsion_strength = 4.0; 
	infection_cell.custom_data.add_variable("signal", "double", 0.0);

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
	scout_cell.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0000002; // 0.2; 
	
	scout_cell.custom_data.add_variable("signal", "double", 0.0);
	scout_cell.custom_data.add_variable("trained", "dimensionless", 0.0);
	scout_cell.custom_data.add_variable("wait_time", "dimensionless", 0.0);

	scout_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.uptake_rates[0] = 0.05; 
	scout_cell.phenotype.secretion.saturation_densities[0] = 1.0; 

	scout_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	scout_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	scout_cell.phenotype.secretion.saturation_densities[0] = 0.0; 


	scout_cell.phenotype.death.rates[apoptosis_model_index] = 0.0 ; 

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
	immune_cell.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.002; // 0.2; 
	
	immune_cell.custom_data.add_variable("wait_time", "dimensionless", 0.0);

	immune_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.uptake_rates[0] = 0.05; 
	immune_cell.phenotype.secretion.saturation_densities[0] = 1.0; 

	immune_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	immune_cell.phenotype.secretion.secretion_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.uptake_rates[0] = 0.0; 
	immune_cell.phenotype.secretion.saturation_densities[0] = 0.0; 

	immune_cell.phenotype.death.rates[apoptosis_model_index] = 0.00001 ; 

	// turn on motility; 
	immune_cell.phenotype.motility.is_motile = true; 
	immune_cell.phenotype.motility.persistence_time = 10.0; 
	immune_cell.phenotype.motility.migration_speed = 0.4;   

	immune_cell.phenotype.mechanics.cell_cell_adhesion_strength = 2.0; 
	immune_cell.phenotype.mechanics.cell_cell_repulsion_strength = 5.0; 

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
	
	std::vector<double> position(3.0, 0.0); 

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
	int signal = 0;

	low = 0.7;
	high = 1.0;

	//Infection cell randomly
	for( int n=0 ; n < number_of_infection ; n++ )
	{
		position[0] = default_microenvironment_options.X_range[0] + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		position[1] = default_microenvironment_options.Y_range[0] + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 

		pC = create_cell(infection_cell); 
		pC->assign_position(position);
		signal = pC->custom_data.find_variable_index("signal");
		pC->custom_data[signal] = (low + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high-low))));
	}

	low = 0.4;
	high = 0.6;

	// Scout cells
	for( int n=0 ; n < number_of_scouts ; n++ )
	{
		position[0] = default_microenvironment_options.X_range[0] + 
				x_range*( relative_outer_margin + (1-2.0*relative_outer_margin)*UniformRandom() ); 
		
		position[1] = default_microenvironment_options.Y_range[0] + 
				y_range*( relative_outer_margin + (1-2.0*relative_outer_margin)*UniformRandom() ); 

		pC = create_cell(scout_cell);
		pC->assign_position(position); 
		signal = pC->custom_data.find_variable_index("signal");
		pC->custom_data[signal] = (low + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high-low))));
	}

	return; 
}


void update_infection_phenotype (Cell* pCell, Phenotype& phenotype, double dt){

		if( phenotype.death.dead == true )
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			return; 
		}

		int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );

		// Sample microenvironment
		static int chem_ind = microenvironment.find_density_index("chemoattractant");
		double chem_conc = pCell->nearest_density_vector()[chem_ind]; 

		static int o2_ind = microenvironment.find_density_index( "oxygen" );
		double o2_conc = pCell->nearest_density_vector()[o2_ind];

		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

		// if (chem_conc<=5.0)
		// {	
		// 	pCell->phenotype.death.rates[apoptosis_model_index] = 0.0001;
		// }

		// else {

			cell_defaults.phenotype.secretion.uptake_rates[o2_ind] = 0.0001;

			// chemoattractant secretion rate
			cell_defaults.phenotype.secretion.secretion_rates[chem_ind] = 50.0;
			cell_defaults.phenotype.secretion.uptake_rates[chem_ind] = 0.0;
			cell_defaults.phenotype.secretion.saturation_densities[chem_ind] = 100.0;

			pCell->phenotype.death.rates[apoptosis_model_index] = 0.00004;
			cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0166; 

			infection_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.1; 
			infection_cell.phenotype.mechanics.cell_cell_repulsion_strength = 20.0;
		//}

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

		static int o2_ind = microenvironment.find_density_index( "oxygen" );
		static int chem_ind = microenvironment.find_density_index( "chemoattractant" );

		double o2_conc = pCell->nearest_density_vector()[o2_ind]; 
		double chem_conc = pCell->nearest_density_vector()[chem_ind]; 

		 if (chem_conc>=10.0) {
			pCell->phenotype.motility.is_motile = false;
		}

		else if(chem_conc<10.0){
			pCell->phenotype.motility.is_motile = true;
			pCell->phenotype.motility.migration_bias = 0.6;
			pCell->phenotype.motility.migration_speed = 0.3;
			phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chem_ind);	
			normalize( &( phenotype.motility.migration_bias_direction ) );
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
	static int o2_ind = microenvironment.find_density_index( "oxygen" );
	double o2_conc = pCell->nearest_density_vector()[o2_ind]; 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	if ((o2_conc < 5.0))
		{
			pCell->phenotype.death.rates[apoptosis_model_index] = 0.00001;
			cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.00001;
		}
}

void update_scout_motility( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true || pCell->is_out_of_domain == true)
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		int chem_index = microenvironment.find_density_index( "chemoattractant" );
		double chem_conc = pCell->nearest_density_vector()[chem_index];
		int mem_index = microenvironment.find_density_index( "mem_attractant" );
		static int t1 = pCell->custom_data.find_variable_index( "trained" );
		static int w = pCell->custom_data.find_variable_index("wait_time");

		if ((chem_conc>0.5) && (pCell->custom_data[t1] == 0.0)){
		pCell->phenotype.motility.migration_bias = 0.5;
		pCell->phenotype.motility.migration_speed = 0.3;
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chem_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );
	}

	return;
}


void scout_rule( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true || pCell->is_out_of_domain == true)
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		Cell* pC = NULL ;
		static int t1 = pCell->custom_data.find_variable_index( "trained" );
		static int w = pCell->custom_data.find_variable_index("wait_time");
		
		int mem_index = microenvironment.find_density_index( "mem_attractant" );
		double time = 0.0;

		std::vector<double> pos = pCell->position;

		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
		
		if(pCell->custom_data[t1] == 0.0)
			{
			for( int n=0 ; n < pCell->cells_in_my_container().size() ; n++ ) {
			pC = pCell->cells_in_my_container()[n]; 

			if( (pC != pCell))
			{
			
			std::vector<double> dist = pC->position;
			dist -= pCell->position;
			double distance = norm(dist);
			int temp = pC->custom_data.find_variable_index("signal");

				if((distance <= 15.0) && (pC->custom_data[temp] >= 0.7))
				{	
				//printf("\nTrained, %d", pCell->ID);
				pCell->custom_data[t1] = 1.0;
				scout_cell.phenotype.mechanics.cell_cell_repulsion_strength = 30.0;
				break;
				} 

			}
		}

		}

		pCell->custom_data[w] = pCell->custom_data[w] + dt;

		if((pCell->custom_data[t1] > 0.8) && (pCell->custom_data[w] >= wait_time)){
			pCell->phenotype.motility.is_motile = true;
			phenotype.motility.migration_bias_direction = pCell->nearest_gradient(mem_index);	
			normalize( &( phenotype.motility.migration_bias_direction ) );
			pCell->phenotype.motility.migration_bias = 0.6;
			pCell->phenotype.motility.migration_speed = 0.8;
			cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.025;
			}

		if ((pos[0] >= x_1-20.0 || pos[0] <= x_0+20.0 || pos[1] >= y_1-20.0 || pos[1] <= y_0+20.0) && (pCell->custom_data[t1] > 0.8)){
				
				#pragma omp critical
				scouts_left += 1;
				printf("\n****** scouts left: %d\n", scouts_left);
				pCell->phenotype.motility.is_motile = false;	
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
	static int o2_ind = microenvironment.find_density_index( "oxygen" );
	double o2_conc = pCell->nearest_density_vector()[o2_ind]; 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	if ((o2_conc < 5.0))
		{
			pCell->phenotype.death.rates[apoptosis_model_index] = 0.00004;
			cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0166;
		}
}

void update_immune_motility( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true || pCell->is_out_of_domain == true)
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		int chem_index = microenvironment.find_density_index( "chemoattractant" );
		double chem_conc = pCell->nearest_density_vector()[chem_index];
		static int w = pCell->custom_data.find_variable_index("wait_time");

		pCell->phenotype.motility.migration_bias = 0.6;
		pCell->phenotype.motility.migration_speed = 0.7;
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chem_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );

	return;
}


void immune_rule( Cell* pCell, Phenotype& phenotype, double dt ){

		if( phenotype.death.dead == true || pCell->is_out_of_domain == true)
		{ 
			pCell->functions.update_migration_bias = NULL; 
			pCell->functions.update_phenotype = NULL;
			pCell->functions.custom_cell_rule = NULL;
			return; 
		}

		Cell* pC = NULL ;
		std::vector<double> pos = pCell->position;

		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

		int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
		
		for( int n=0 ; n < pCell->cells_in_my_container().size() ; n++ ) {
		pC = pCell->cells_in_my_container()[n]; 

		if( (pC != pCell))
		{
		
		std::vector<double> dist = pC->position;
		dist -= pCell->position;
		double distance = norm(dist);
		int temp = pC->custom_data.find_variable_index("signal");

			if((distance <= 15.0) && (pC->custom_data[temp] >= 0.7))
			{	
			pC->start_death(apoptosis_model_index);
			pCell->phenotype.death.rates[apoptosis_model_index] = 0.0001;
			} 

		}
	}

	return;
}


// Defining my own custom coloring function
std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// color based on alive and dead cells but can start 
	// with the simple coloring scheme. 
	
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
		output[1] = "cyan";
		output[3] = "cyan";
		output[0] = "cyan"; 
		output[2] = "cyan";

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

	static int chem_index = microenvironment.find_density_index( "chemoattractant" );
	double chem_conc = 0.0;

	int scout_cells, new_untrained_cells, new_trained_cells = 0;
	int number_of_scouts = parameters.ints("number_of_scouts"); 

	// collect the chemoattractant concentration from the domain (all voxels)
	#pragma omp parallel for 
	for( int k=0 ; k < microenvironment.number_of_voxels(); k++ )
	{
	chem_conc += microenvironment(k)[chem_index];
	}

	//printf("time stamp:,%f", dt);

	//# pragma omp critical
	if ( (chem_conc >= 200.0) & (scouts_left > 0) ) {

		printf("\nConc. %f , sending in the troops.\n",chem_conc);
		
		//spawn cells
		scout_cells = scouts_left + (dt *  ( (chem_conc  * 0.4 * scouts_left) - (0.05 * scouts_left) - (0.0015 * scouts_left * untrained_immune_cell) ) );
		
		new_untrained_cells = new_untrained_cells +  (dt * ( (chem_conc * 0.04 * untrained_immune_cell) - (0.0005 * untrained_immune_cell) - ( 0.00015 * scout_cells * untrained_immune_cell) ) );

		new_trained_cells = new_trained_cells + ( dt * ( (chem_conc * 0.004 * trained_immune_cell) - (0.005 * trained_immune_cell) + (0.15  * scout_cells * untrained_immune_cell) ) - (0.00015 * trained_immune_cell) );
		
		printf("\n scout cells: %d, untrained cells: %d, trained cells %d\n", scout_cells, new_untrained_cells, new_trained_cells);
		// call the setup tissue function again.
		untrained_immune_cell = new_untrained_cells;
		trained_immune_cell = new_trained_cells;
		scouts_left = 0;
		setup_new_cells(scout_cells, new_trained_cells);

	}
	return;
}

void setup_new_cells(int scout_cells, int trained_cells){

	int number_of_trained = trained_cells;   
	int number_of_scouts = scout_cells;  

	std::cout << "Placing new cells ... " << std::endl; 
	
	// randomly place seed cells 
	
	std::vector<double> position(5.0, 0.0); 

	x_0 = default_microenvironment_options.X_range[0]; 
	x_1 = default_microenvironment_options.X_range[1];

	y_0 = default_microenvironment_options.Y_range[0]; 
	y_1 = default_microenvironment_options.Y_range[1];

	double x_range = x_1 - x_0; 
	double y_range = y_1 - y_0; 

	double relative_margin = 0.5;  
	double relative_outer_margin = 0.04;

	double xperiod = 30.0;
	double yperiod = 20.0;

	low = 0.5;
	high = 0.6;
	int signal = 0;

	Cell* pC;

	#pragma omp critical
	for( int n=0 ; n < number_of_scouts ; n++ )
	{
		// position[0] = default_microenvironment_options.X_range[0] + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		// position[1] = default_microenvironment_options.Y_range[0] + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 

		//#pragma omp critical
		position[0] = default_microenvironment_options.X_range[0] + xperiod;
		position[1] = default_microenvironment_options.Y_range[0] + yperiod;
		pC = create_cell(scout_cell);
		pC->assign_position(position); 
		signal = pC->custom_data.find_variable_index("signal");
		pC->custom_data[signal] = (low + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high-low))));
		xperiod += xperiod;
	}

	xperiod = 20.0;

	for( int n=0 ; n < number_of_trained ; n++ )
	{
		// position[0] = default_microenvironment_options.X_range[0] + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		// position[1] = default_microenvironment_options.Y_range[0] + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() );  
		position[0] = default_microenvironment_options.X_range[0] + xperiod;
		position[1] = default_microenvironment_options.Y_range[1] - yperiod;
		pC = create_cell(immune_cell);
		pC->assign_position(position); 
		xperiod += xperiod;
	}

	return;
}
