#include "MPIBodyContainer.hpp"
YADE_PLUGIN((MPIBodyContainer));
CREATE_LOGGER(MPIBodyContainer);


void MPIBodyContainer::clearContainer() {
      bContainer.clear(); 
}
