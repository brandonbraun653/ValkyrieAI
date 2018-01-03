#pragma once
#ifndef GA_CONFIG_H_
#define GA_CONFIG_H_

/*-----------------------------------------------
* Global Genetic Algorithm "Stuff"
*-----------------------------------------------*/
#define DEBUGGING_ENABLED					0 
#define MAX_THREADS_ALLOWED_PER_OPTIMIZER	16

/*-----------------------------------------------
* Runtime Options
*-----------------------------------------------*/
/* Initialize Population Function */
#define GA_TRACE_INITIALIZE_POPULATION		0

/* Evaluate Model Function */
#define GA_TRACE_EVALUATE_MODEL				0 
#define FIXED_PID_VALUES					0 

/* Population Filter Function */
#define GA_TRACE_FILTER_POPULATION			0 

/* Evaluate Fitness Function */
#define GA_TRACE_EVALUATE_FITNESS			0 
#define GA_TRACE_SELECT_PARENTS				0 

/* Breed Generation Function */
#define GA_TRACE_BREED_GENERATION			0 
#define GA_ENFORCE_RESOLUTION_BG			0 

/* Mutate Generation Function */	
#define GA_TRACE_MUTATE_GENERATION			0 
#define GA_ENFORCE_RESOLUTION_MG			0 

/* Check Convergence Function */
#define GA_TRACE_CHECK_CONVERGENCE			0 
#define GA_REPORT_DATA_CHECK_CONVERGENCE	0 


/* Debugging Statements */
#define SS_TRACE							0
#define SS_TRACE_LOG_COUT					0
#define SS_TRACE_LOG_CSV					0
#define SS_THREAD_TRACE_LOG_COUT			0
#define SS_THREAD_TRACE_LOG_CSV				0


#endif