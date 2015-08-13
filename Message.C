#include "Message.H"



Message::Message(int numParms) 
  : parms(numParms), numParms(numParms)
{
  computeBufferSize();
  buffer=(void*) new char[bufferSize];
  firstSlot=buffer+sizeof(MessageType);
}



void Message::computeBufferSize()
{
  bufferSize=
    sizeof(MessageType)       // tag
    +sizeof(double)*numParms  // parms
    +2*sizeof(long)           // begin & end
    +sizeof(int);             // order
}



Message::~Message() 
{
  delete [] buffer;
}



void Message::resize(int numParms) {
  parms=numParms;
  this->numParms=numParms;
  delete [] buffer;
  computeBufferSize();
  buffer=(void*) new char[bufferSize];
  firstSlot=buffer+sizeof(MessageType);
}



void Message::packLikelihood(double likelihood)
{
  packType(MSG_LIKELIHOOD);
  memcpy(firstSlot,(void*)&likelihood,sizeof(likelihood));
}



void Message::packTreeNotFixed() {
  packType(MSG_TREE_NOT_FIXED);
}



double Message::getLikelihood() 
{
  return likelihood;
}



int Message::getBufferSize()
{
  return bufferSize;
}



void *Message::getBuffer() 
{
  return buffer;
}



void Message::packInterval(long begin,long end) 
{
  packType(MSG_SET_INTERVAL);
  void *slot=firstSlot;
  memcpy(slot,(void*)&begin,sizeof(begin));
  slot+=sizeof(begin);
  memcpy(slot,(void*)&end,sizeof(end));
}



void Message::packEvaluate(const GSL::Vector &parms) 
{
  packType(MSG_EVALUATE);
  void *slot=firstSlot;
  for(int i=0 ; i<numParms ; ++i) {
    double elem=parms[i];
    memcpy(slot,(void*)&elem,sizeof(elem));
    slot+=sizeof(elem);
  }
}



void Message::packQuit() 
{
  packType(MSG_QUIT);
}



void Message::packType(MessageType type) 
{
  memcpy(buffer,(void*)&type,sizeof(type));
}



MessageType Message::unpack() 
{
  void *slot;
  double parm;
  MessageType type;
  memcpy((void*)&type,buffer,sizeof(type));
  switch(type) 
    {
    case MSG_SET_INTERVAL:
      slot=firstSlot;
      memcpy((void*)&begin,slot,sizeof(begin));
      slot+=sizeof(begin);
      memcpy((void*)&end,slot,sizeof(end));
      break;
    case MSG_EVALUATE:
      slot=firstSlot;
      for(int i=0 ; i<numParms ; ++i) {
	memcpy((void*)&parm,slot,sizeof(parm));
	parms[i]=parm;
	slot+=sizeof(parm);
      }
      break;
    case MSG_LIKELIHOOD:
      memcpy((void*)&likelihood,firstSlot,sizeof(likelihood));
      break;
    case MSG_NEXT_HIGHER_ORDER:
      memcpy((void*)&order,firstSlot,sizeof(order));
      slot=firstSlot+sizeof(order);
      for(int i=0 ; i<numParms ; ++i) {
	memcpy((void*)&parm,slot,sizeof(parm));
	parms[i]=parm;
	slot+=sizeof(parm);
      }
      break;
    case MSG_TREE_NOT_FIXED:
      break;
    case MSG_QUIT: 
      break;
    }
  return type;
}



GSL::Vector &Message::getParms() 
{ 
  return parms; 
}



long Message::getBegin() 
{ 
  return begin; 
}



long Message::getEnd() 
{ 
  return end; 
}



void Message::packNextHigherOrder(int order,const GSL::Vector &parms) {
  packType(MSG_NEXT_HIGHER_ORDER);
  memcpy(firstSlot,(void*)&order,sizeof(order));

  void *slot=firstSlot+sizeof(order);
  for(int i=0 ; i<numParms ; ++i) {
    double elem=parms[i];
    memcpy(slot,(void*)&elem,sizeof(elem));
    slot+=sizeof(elem);
  }
}



int Message::getOrder() {
  return order;
}
