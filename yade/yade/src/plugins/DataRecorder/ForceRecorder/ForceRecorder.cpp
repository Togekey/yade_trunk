#include "ForceRecorder.hpp"
#include "RigidBodyParameters.hpp"
#include "Omega.hpp"
#include "ComplexBody.hpp"
#include "ActionParameterForce.hpp"

#include <boost/lexical_cast.hpp>

ForceRecorder::ForceRecorder () : Actor(), actionForce(new ActionParameterForce)
{
	outputFile = "";
	interval = 50;
	startId = 0;
	endId = 1;
	changed = false;
}

void ForceRecorder::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
		ofile.open(outputFile.c_str());
	}
}

void ForceRecorder::registerAttributes()
{
	Actor::registerAttributes();
	REGISTER_ATTRIBUTE(outputFile);
	REGISTER_ATTRIBUTE(interval);
	REGISTER_ATTRIBUTE(startId);
	REGISTER_ATTRIBUTE(endId);
	REGISTER_ATTRIBUTE(bigBallId);
	REGISTER_ATTRIBUTE(bigBallReleaseTime);
}

bool ForceRecorder::isActivated()
{
	return ((Omega::instance().getCurrentIteration() % interval == 0) && (ofile));
}

void ForceRecorder::action(Body * body)
{
	ComplexBody * ncb = dynamic_cast<ComplexBody*>(body);
	Real x=0, y=0, z=0;
	
	for( unsigned int i = startId ; i <= endId ; ++i )
	{
		Vector3r force = dynamic_cast<ActionParameterForce*>(ncb->actions->find( i , actionForce->getClassIndex() ) . get() )->force;
		
		x+=force[0];
		y+=force[1];
		z+=force[2];
	}
	
	ofile << lexical_cast<string>(Omega::instance().getSimulationTime()) << " " 
		<< lexical_cast<string>(x) << " " 
		<< lexical_cast<string>(y) << " " 
		<< lexical_cast<string>(z) << endl;

		
	// FIXME all that lines do not belong to ForceRecorder
	if( bigBallReleaseTime < Omega::instance().getSimulationTime() && (!changed) )
	{
		changed = true;
		(*(ncb->bodies))[bigBallId]->isDynamic = true;
	}
}

