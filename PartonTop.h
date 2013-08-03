#ifndef PARTONTOP__HH 1
#define PARTONTOP__HH 1

#include "fjClustering.h"
#include "MyTopEvent.h"


class PartonTop
{
public:

    //constuctor. send in  the ttbar pair
	PartonTop(objets &ttbar);

	DecayChannels( pair<int,int> &);
  


private:

	Pythia8::Particle m_top, m_tbar;


};

#endif
