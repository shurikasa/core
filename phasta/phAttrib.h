#ifndef PH_ATTRIB_H
#define PH_ATTRIB_H

#include "phBC.h"

namespace ph {

void getSimmetrixAttributes(gmi_model* model, BCs& bcs);
void clearAttAssociation(gmi_model* mdl, ph::Input& in);

}

#endif
