#include "valkyrie_engine.h"

using namespace boost::interprocess;

/*-----------------------------------------------
* Constructor/Destructor
*-----------------------------------------------*/
ValkyrieEngine::ValkyrieEngine()
{
	ostream_mtx_sPtr = boost::make_shared<boost::mutex>();

}

ValkyrieEngine::~ValkyrieEngine()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
FCSOptimizer_Handle ValkyrieEngine::newOptimizer(FCSOptimizer_Init_t init)
{
	FCSOptimizer_Handle optimizer = boost::make_shared<FCSOptimizer_Handle_t>();

	optimizer->Init = init;
	optimizer->Engine = boost::make_shared<FCSOptimizer>();
	optimizer->Status = GA_SETUP;

	/* Keep a reference to the new instance */
	optimizerInstances.push_back(optimizer);

	return optimizer;
}

void ValkyrieEngine::initialize(const FCSOptimizer_Handle& optimizer)
{
	optimizer->Engine->init(optimizer->Init);

	
	/* The Engine's command queue will have been initialized now, so setup
	the interface queue for the main thread. */
	try
	{
		optimizer->CommandQueue = new message_queue(
			open_only,
			optimizer->Init.messageQueueName.data()
		);
	}
	catch (interprocess_exception &ex)
	{
		//TODO: Use the console mutex here to make sure the output is readable 
		std::cout << ex.what() << " in the main thread." << std::endl;
		std::cout << "Cannot properly open the " 
			<< optimizer->Init.optimizerName 
			<< " Optimizer message queue. Unable to send commands." 
			<< std::endl;
	}

	/* Start the thread, but it won't execute the optimization process until 
	the start command is given from the user. */
	auto threadFunc = boost::bind(&FCSOptimizer::run, optimizer->Engine);

	optimizer->Thread = boost::make_shared<boost::thread>(threadFunc);
}

void ValkyrieEngine::start(const FCSOptimizer_Handle& optimizer)
{
	#if defined(_WIN32) && !defined(_WIN64)
	/* Start up the TCP connection before continuing */
	if (optimizer->Init.solverParam.modelType == GA_MODEL_NEURAL_NETWORK)
	{
		NN_TCPModel_sPtr nn_model = boost::dynamic_pointer_cast<NN_TCPModel, NN_ModelBase>(optimizer->Init.neuralNetModel);
		nn_model->openConnection();
	}
	#endif

	int command = START;
	optimizer->CommandQueue->try_send(&command, sizeof(command), 0);
}

void ValkyrieEngine::pause(const FCSOptimizer_Handle& optimizer)
{
	int command = PAUSE;
	optimizer->CommandQueue->try_send(&command, sizeof(command), 0);
}

void ValkyrieEngine::stop(const FCSOptimizer_Handle& optimizer)
{
	int command = STOP;
	optimizer->CommandQueue->try_send(&command, sizeof(command), 0);
}

void ValkyrieEngine::initAll()
{
	for (int i = 0; i < optimizerInstances.size(); i++)
		initialize(optimizerInstances[i]);
}

void ValkyrieEngine::startAll()
{
	for (int i = 0; i < optimizerInstances.size(); i++)
		start(optimizerInstances[i]);
}

void ValkyrieEngine::pauseAll()
{
	for (int i = 0; i < optimizerInstances.size(); i++)
		pause(optimizerInstances[i]);
}

void ValkyrieEngine::stopAll()
{
	for (int i = 0; i < optimizerInstances.size(); i++)
		stop(optimizerInstances[i]);
}

void ValkyrieEngine::waitForCompletion()
{
	/* Assuming a thread exists and is joinable, force joining before exit. */
	for (int i = 0; i < optimizerInstances.size(); i++)
	{
		if (optimizerInstances[i]->Thread)
		{
			if (optimizerInstances[i]->Thread->joinable())
				optimizerInstances[i]->Thread->join();
		}
	}
}
